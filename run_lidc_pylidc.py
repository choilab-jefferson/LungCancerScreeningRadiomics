#!/usr/bin/env python3
import os.path as osp
from multiprocessing import freeze_support

import ruffus
import ruffus.cmdline as cmdline

import utils
import tasks

##########################################################################
# Set your environmental parameters
##########################################################################
# pylidc - https://pylidc.github.io/install.html

output_dir = 'output'
data_dir = 'DATA'
temp_dir = 'tmp'
iso_size = '1'

##########################################################################


def make_pipeline_LIDC_analysis(experiment_set):
    pipeline_name = experiment_set
    data_path = osp.join(data_dir, experiment_set)
    output_path = osp.join(output_dir, experiment_set)
    feature_list_path = osp.join(
        output_dir, "feature-list_" + experiment_set + ".csv")
    pipeline = ruffus.Pipeline(pipeline_name)
    
    starting_file_names = tasks.load_scan_list(data_path)
    #print(starting_file_names)

    pipeline.originate(name="task_originate",
                       task_func=tasks.originate_pylidc,
                       output=starting_file_names) \
        .follows(ruffus.mkdir(output_path)) \
        .follows(ruffus.mkdir(data_path)) \
        .follows(ruffus.mkdir(data_path+"/CT"))

    pipeline.transform(name="task_image_resample",
                       task_func=tasks.image_resample,
                       input=ruffus.output_from("task_originate"),
                       filter=ruffus.formatter(),
                       output="{path[0]}/{basename[0]}-" + iso_size + "mm.nrrd")

    pipeline.subdivide(name="task_check_nodules",
                       task_func=tasks.check_nodules,
                       input=ruffus.output_from("task_originate"),
                       filter=ruffus.formatter(
                           experiment_set + r"/CT/(?P<pid>[^/]+)_CT.nrrd"),
                       output=["{subpath[0][1]}/{pid[0]}/{basename[0]}_*.nodule"])

    pipeline.transform(name="task_extract_nodule_labels",
                       task_func=tasks.extract_nodule_labels,
                       input=ruffus.output_from("task_check_nodules"),
                       filter=ruffus.formatter(),
                       output=["{path[0]}/{basename[0]}-Phy1-label.nrrd",
                               "{path[0]}/{basename[0]}-Phy2-label.nrrd",
                               "{path[0]}/{basename[0]}-Phy3-label.nrrd",
                               "{path[0]}/{basename[0]}-Phy4-label.nrrd"])

    pipeline.subdivide(name="task_feature_extraction_phy",
                       task_func=tasks.feature_extraction_phy,
                       input=ruffus.output_from("task_extract_nodule_labels"),
                       filter=ruffus.formatter(
                           experiment_set + r"/(?P<pid>[^/]+)/(?P<image_name>.*)_(?P<nid>\d*)-.*-label.nrrd"),
                       output=[output_dir + "/" + experiment_set + "/{image_name[0]}_{nid[0]}-Phy1-" + iso_size + "mm.txt",
                               output_dir + "/" + experiment_set + "/{image_name[0]}_{nid[0]}-Phy2-" + iso_size + "mm.txt",
                               output_dir + "/" + experiment_set + "/{image_name[0]}_{nid[0]}-Phy3-" + iso_size + "mm.txt",
                               output_dir + "/" + experiment_set + "/{image_name[0]}_{nid[0]}-Phy4-" + iso_size + "mm.txt"],
                       extras=["{subpath[0][1]}/CT/{image_name[0]}.nrrd"])\
        .follows("task_image_resample")

    pipeline.transform(name="task_staple_comparison",
                       task_func=tasks.staple_comparison,
                       input=ruffus.output_from("task_extract_nodule_labels"),
                       filter=ruffus.formatter(
                           experiment_set + r"/(?P<pid>[^/]+)/(?P<image_name>.*)_(?P<nid>\d*)-.*-label.nrrd"),
                       output=["{path[0]}/{image_name[0]}_{nid[0]}-all-label.nrrd",
                               "{path[0]}/{image_name[0]}_{nid[0]}-STAPLE-label.nrrd",
                               output_dir + "/" + experiment_set + "/{image_name[0]}_{nid[0]}-STAPLE.txt"])

    pipeline.transform(name="task_segment_nodule",
                       task_func=tasks.segment_nodule,
                       input=ruffus.output_from("task_staple_comparison"),
                       filter=ruffus.formatter(
                           experiment_set + r"/(?P<pid>[^/]+)/(?P<image_name>.*)_(?P<nid>\d*)-all-label.nrrd"),
                       output=["{path[0]}/{image_name[0]}_{nid[0]}-seg-label.nrrd",
                               "{path[0]}/{image_name[0]}_{nid[0]}-levelset.nrrd"],
                       extras=["{subpath[0][1]}/CT/{image_name[0]}.nrrd"])\

    pipeline.transform(name="task_feature_extraction_seg",
                       task_func=tasks.feature_extraction,
                       input=ruffus.output_from("task_segment_nodule"),
                       filter=ruffus.formatter(
                           experiment_set + r"/(?P<pid>[^/]+)/(?P<image_name>.*)_(?P<nid>\d*)-seg-label.nrrd"),
                       output=output_dir + "/" + experiment_set +
                       "/{image_name[0]}_{nid[0]}-seg-" + iso_size + "mm.txt",
                       extras=["{subpath[0][1]}/CT/{image_name[0]}.nrrd"])\
        .follows("task_image_resample")

    pipeline.merge(name="task_staple_organization",
                   task_func=tasks.staple_organization,
                   input=ruffus.output_from("task_staple_comparison"),
                   output=output_dir + "/STAPLE_" + experiment_set + ".csv")

    pipeline.transform(name="task_feature_extraction",
                       task_func=tasks.feature_extraction,
                       input=ruffus.output_from("task_staple_comparison"),
                       filter=ruffus.formatter(
                           experiment_set + r"/(?P<pid>[^/]+)/(?P<image_name>.*)_(?P<nid>\d*)-all-label.nrrd"),
                       output=output_dir + "/" + experiment_set +
                       "/{image_name[0]}_{nid[0]}-all-" + iso_size + "mm.txt",
                       extras=["{subpath[0][1]}/CT/{image_name[0]}.nrrd"])\
        .follows("task_image_resample")

    pipeline.merge(name="task_feature_organization",
                   task_func=tasks.feature_organization,
                   input=ruffus.output_from(
                       ["task_feature_extraction", "task_feature_extraction_phy", "task_feature_extraction_seg"]),
                   output=feature_list_path)

    return pipeline


if __name__ == "__main__":
    freeze_support()

    experiment_set = 'nodule-lidc_pylidc'
    pipeline_LIDC = make_pipeline_LIDC_analysis(experiment_set)

    parser = cmdline.get_argparse(description='LIDC radiomics')

    options = parser.parse_args()

    # standard python logger which can be synchronised across concurrent
    # Ruffus tasks
    logger, logger_mutex = cmdline.setup_logging(
        __name__, options.log_file, options.verbose)

    #
    #   Run
    #
    cmdline.run(options)
