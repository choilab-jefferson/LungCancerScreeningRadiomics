%% main function
clc;
clear;
close all

set_environment; % call environment setting

%% spiculation toolbox
addpath(genpath('util'))


%% set values
experiment_set = 'MSKLung';

%% directory paths
experiment_path = [data_path '/' experiment_set];
smooth = true;
if smooth
    obj_path=[experiment_path '/objs_smooth'];
    sph_map_path = [experiment_path '/spherical_obj_smooth'];
else
    obj_path=[experiment_path '/objs'];
    sph_map_path = [experiment_path '/spherical_obj'];
end

if ~isdir(sph_map_path); mkdir(sph_map_path); end
dir_obj = dir([obj_path '/*.obj']);
pid_list = {dir_obj.name};

features = table();
all_spikes = table();
merged_nodule_info = [];

%% main process
parfor idx = 1:size(pid_list,2)
    tic
    try
        id = strsplit(pid_list{idx},'.');
        id = id{1};
        ids = strsplit(id,'_');
        pid = ids{1}; nid = ids{2};
        nids = strsplit(nid,'-'); % to cover multiple segmentations for a nodule

        fprintf('%d %s %s\n', idx, pid, nid);

        %if numel(strfind(pid, '_'))>0 || numel(strfind(pid, 'Ell'))>0 || ~strcmp(nid,'1'), continue; end
        %if numel(strfind(pid, '0043'))==0, continue; end


        %%
        obj_filename = [obj_path '/' pid '_' nid '.obj'];
        sph_map_filename = [sph_map_path '/' pid '_' nid '_spherical.obj'];

        spherical_param(sph_map_filename, obj_filename)
    catch
        continue
    end
    toc
end