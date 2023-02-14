%% main function
clc;
clear;

set_environment; % call environment setting


%% set values
sz_bbox = 200;
smooth = true;
experiment_set = 'MSKLung';


%% directory paths
experiment_path=[data_path '/' experiment_set];
if smooth
    objs_experiment_path=[experiment_path '/objs_smooth'];
else
    objs_experiment_path=[experiment_path '/objs'];
end
    
if ~isdir(objs_experiment_path); mkdir(objs_experiment_path); end

dir_nrrd = dir(experiment_path);
pid_list = {dir_nrrd.name};

features = table();
merged_nodule_info = table();
%% main process
for idx = 1:numel(pid_list)
    tic
    pid = pid_list{idx};

    if sum(strfind(pid,'.') < 3) || sum(strfind(pid,'objs')) > 0 || strcmp(pid,'CT')
        continue
    end
    %idif numel(strfind(pid, '0020'))==0, continue; end

%     if sum(strcmp(selected, pid)) == 0
%         continuecp -LR results/MSKLung /data/apps/users/wxc151/MSKLung

%     end


    results = table();

    fprintf('%d %s\n', idx, pid);
    label_files = dir([experiment_path '/' pid '/' pid '_CT_*-R1' iso_t '-label.nrrd']);
    for l = 1:numel(label_files)
        tic
        label_file = label_files(l);
        %% input part
        label_file_path = [label_file.folder '/' label_file.name];
        id = strsplit(label_file.name,'_');
        pid = id{1};
        id = strsplit(id{3},'-');
        nid = [id{1} '-' id{2}];
        model_file_path = [objs_experiment_path '/' pid '_' nid '.obj'];
        try
            s = generate_mesh_model(model_file_path, label_file_path, smooth);
        catch
            continue
        end
        toc
    end
end

%writetable(merged_nodule_info, 'nodule_info.csv')

