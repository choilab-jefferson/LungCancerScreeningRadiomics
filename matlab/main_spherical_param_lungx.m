%% main function
clc;
clear;
close all

set_environment; % call environment setting

%% spiculation toolbox
addpath(genpath('util'))


%% set values
experiment_set = 'nodule-lungx';

%% directory paths
experiment_path = [data_path '/' experiment_set];
obj_path = [experiment_path '/objs'];
sph_map_path = [experiment_path '/spherical_obj'];

if ~isdir(sph_map_path); mkdir(sph_map_path); end
dir_obj = dir([obj_path '/*.obj']);
pid_list = {dir_obj.name};

features = table();
all_spikes = table();
merged_nodule_info = [];
%% main process
for idx = 1:size(pid_list,2)
    tic
    id = strsplit(pid_list{idx},'.');
    id = id{1};
    ids = strsplit(id,'_');
    pid = ids{1}; nid = ids{2};
    nids = strsplit(nid,'-'); % to cover multiple segmentations for a nodule
    
    fprintf('%d %s %s\n', idx, pid, nid);
    
    %if numel(strfind(pid, '_'))>0 || numel(strfind(pid, 'Ell'))>0 || ~strcmp(nid,'1'), continue; end
    %if numel(strfind(pid, '0043'))==0, continue; end
    
    %%
    if ~exist('refine','var')
        refine = false;
    end
    if ~exist('smooth','var')
        smooth = false;
    end
    
    %%
    if ~isempty(nid)
        obj_filename = [obj_path '/' pid '_' nid '.obj'];
        sph_map_filename = [sph_map_path '/' pid '_' nid '_spherical.obj'];
    else
        obj_filename = [obj_path '/' pid '.obj'];
        sph_map_filename = [sph_map_path '/' pid '_spherical.obj'];
    end
    
    spherical_param(sph_map_filename, obj_filename)
    toc
end