%% main function
clc;
clear;

set_environment; % call environment setting


%% set values
% selected = selected_72;
% experiment_set = 'nodule-radiomics_72';

experiment_set = 'nodule-lidc';
sz_bbox = 200;
smooth = true;


%% directory paths
experiment_path=[data_path '/' experiment_set];
objs_experiment_path=[experiment_path '/objs'];
if ~isdir(objs_experiment_path); mkdir(objs_experiment_path); end

dir_nrrd = dir(experiment_path);
pid_list = {dir_nrrd.name};

features = table();
merged_nodule_info = table();
%% main process
for idx = 1:numel(pid_list)
    tic
    pid = pid_list{idx};

    if numel(strfind(pid, 'LIDC')) == 0
        continue
    end
    
%     if numel(strfind(pid, '0020'))==0, continue; end

%     if sum(strcmp(selected, pid)) == 0
%         continue
%     end


    results = table();

    tic % tic starts a stopwatch timer
    fprintf('%d %s\n', idx, pid);

    nodule_info = readtable([experiment_path '/' pid '/' pid '_all.csv']);
    nodule_info(isnan(nodule_info.Volume),:) =[];
    nodule_info = sortrows(nodule_info, 'Volume', 'descend');

    if ~isdir(['output/' pid]); mkdir(['output/' pid]); end
    for nid = 1:size(nodule_info,1)
        rnid = num2str(nodule_info{nid, 'nid'});
        label_file_path = [experiment_path '/' pid '/' pid '_CT_' rnid '-seg' iso_t '-label.nrrd'];
        model_file_path = [objs_experiment_path '/' pid '_' rnid '.obj'];
        try
            s = generate_mesh_model(model_file_path, label_file_path, smooth);
        catch exception
            disp(exception)
        end
        toc
    end
end

writetable(merged_nodule_info, [objs_experiment_path '/nodule_info.csv'])

