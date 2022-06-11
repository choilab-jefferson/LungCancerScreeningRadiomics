%% main function
clc;
clear;
close all

set_environment; % call environment setting

%% spiculation toolbox
addpath(genpath('util'))


%% set values
experiment_set = 'nodule-lungx';
esph_factor = (3/(4*pi))^(1/3);


%% directory paths
subset = '';
experiment_path = [data_path '/' experiment_set];

if strfind(subset,'seg')
    obj_path = [experiment_path '/objs_seg'];
    sph_map_path = [experiment_path '/spherical_obj_seg'];
else
    obj_path = [experiment_path '/objs'];
    sph_map_path = [experiment_path '/spherical_obj'];
end
if ~isdir(sph_map_path); mkdir(sph_map_path); end
dir_obj = dir([obj_path '/*.obj']);
pid_list = {dir_obj.name};

features = table();
all_peaks = table();
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
    [s.vertices,s.faces] = readOBJ(obj_filename);
    if refine, s = refinepatch(s); end
    if smooth, s = smoothpatch(s,1,1); end
    temp_ply = [tempdir experiment_set subset '.ply'];
    writePLY(temp_ply, s.vertices, s.faces);
    disp('Spherical parameterization...')
    system([ConformalizedMCF ' --in ' temp_ply ' --outHeader ' tempdir ' --steps 3000 --threads 1']);
    [s1.vertices,s1.faces] = readPLY([tempdir '.3000.ply']);
    writeOBJ(sph_map_filename, s1.vertices, s1.faces);
    toc
end