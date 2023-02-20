%% main function
clc;
clear;
close all

set_environment; % call environment setting

%% spiculation toolbox
addpath(genpath('utils'))
addpath(genpath('spiculation_detection'))

%% set values
experiment_set = 'MSKLung';


data_path = '/home/wxc151/gitRepos/radiomics_pipelines/MSKLung_pipeline/results';

%% directory paths
experiment_path = [data_path '/' experiment_set];
smooth = false;
if smooth
    obj_path=[experiment_path '/objs_smooth'];
    sph_map_path = [experiment_path '/spherical_obj_smooth'];
    output_experiment_path=['../' output_path '/' experiment_set '_smooth'];
    smooth_t = '-smooth';
else
    obj_path=[experiment_path '/objs'];
    sph_map_path = [experiment_path '/spherical_obj'];
    output_experiment_path=['../' output_path '/' experiment_set];
    smooth_t = '';
end

if ~isdir(output_experiment_path); mkdir(output_experiment_path); end
if ~isdir([output_experiment_path '/spiculation']); mkdir([output_experiment_path '/spiculation']); end
if ~isdir([output_experiment_path '/parameters']); mkdir([output_experiment_path '/parameters']); end


obj_path = experiment_path;
sph_map_path = experiment_path;

dir_obj = dir([obj_path '/*/*_spherical.obj']);
pid_list = {dir_obj.name};

features = table();
all_spikes = table();
merged_nodule_info = [];
n_info = [];

%% main process
for idx = 1:size(pid_list,2)
    tic % tic starts a stopwatch timer

    %% extract PID and NID from file name
    id = strsplit(pid_list{idx},'.');
    filename = strrep(id{1}, '_spherical', '');
    id = id{1};
    ids = strsplit(id,'_');
    pid = ids{1}; nid = ids{3};
    nids = strsplit(nid,'-');
    nid = [nids{2} '-' nids{3}];
    iso_t = ['-' strrep(nids{4}, 'seg_', '')];
    
%     if ~strcmp(pid, '190365') || ~strcmp(nid, '1-1')
%         continue
%     end

    fprintf('%d %s %s\n', idx, pid, nid);

    %% ignore some cases
    if numel(strfind(pid, '_'))>0 || numel(strfind(pid, 'Ell'))>0, continue; end
    %if numel(strfind(pid, '_'))>0 || numel(strfind(pid, 'Ell'))>0 || ~strcmp(nid,'1'), continue; end    
    %if numel(strfind(pid, '0289'))==0, continue; end


    %% input & ouput files
    label_file = [experiment_path '/' pid '/' strrep(filename, '-py', '') '-label.nrrd'];
    attached_label_filename = [experiment_path '/' pid '/' pid '_CT_R1-' nid iso_t '-attached-label.nrrd'];
    obj_filename = [experiment_path '/' pid '/' filename '.obj'];
    sph_map_filename = [experiment_path '/' pid '/' filename '_spherical.obj'];
    
    ard_voxel_filename = [output_experiment_path '/'  pid '/' pid '_CT_' nid iso_t '-R1-ard.nrrd'];
    ard_surface_filename = [output_experiment_path '/'  pid '/' pid '_CT_' nid iso_t '-R1-ard-surface.nrrd'];
    spike_lable_filename = [output_experiment_path '/'  pid '/' pid '_CT_' nid iso_t '-R1-spikes-label.nrrd'];
    spike_surface_filename = [output_experiment_path '/'  pid '/' pid '_CT_' nid iso_t '-R1-spikes-surface-label.nrrd'];


    %% additional output
    ard_filename = [output_experiment_path '/parameters/' filename '_nd.txt'];
    attached_indices_filename = [output_experiment_path '/parameters/' filename '_attached.txt'];
    spike_apex_filename = [output_experiment_path '/parameters/' filename '_apex.txt'];
    spiculation_fcsv_filename = [output_experiment_path '/parameters/' filename '.fcsv'];
    feature_filename = [output_experiment_path '/parameters/' filename '_features.txt'];
    figure_fig_filename = [output_experiment_path '/spiculation/' filename '_sphere_param.fig'];
    figure_png_filename = [output_experiment_path '/spiculation/' filename '_sphere_param.png'];

    if ~isdir([output_experiment_path '/' pid]); mkdir([output_experiment_path '/' pid]); end
   
    try
        [f, spikes_table]  = spiculation_pipeline(pid, nid, filename, n_info, ...
            label_file, attached_label_filename, obj_filename, sph_map_filename, ...
            ard_voxel_filename, ard_surface_filename, spike_lable_filename, spike_surface_filename, ...
            ard_filename, attached_indices_filename, spike_apex_filename, spiculation_fcsv_filename, feature_filename, ...
            figure_fig_filename, figure_png_filename);
        all_spikes = [all_spikes; spikes_table];
        features = [features;f];
    catch exception
        disp(exception)
        %continue
    end

    clear colors
    close all
    toc
end
writetable(all_spikes, [output_experiment_path '/spiculation/spikes.csv']);
writetable(features, [output_experiment_path '/spiculation/features.csv']);
