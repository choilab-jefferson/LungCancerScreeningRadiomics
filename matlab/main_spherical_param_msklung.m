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

        s = struct();
        s1 = s;
        [s.vertices,s.faces] = readOBJ(obj_filename);
        tempheader = tempname;
        mkdir(tempheader)
        temp_ply = [tempheader '/' experiment_set '_' id '_' nid '.ply'];
        out_header = [tempheader '/out'];
        writePLY(temp_ply, s.vertices, s.faces);
        disp('Spherical parameterization...')
        ConformalizedMCF = ['docker run --user 1007 -v ' tempheader ':' tempheader ' wookjinchoi/conformalized_mcf:latest ConformalizedMCF'];
        system([ConformalizedMCF ' --in ' temp_ply ' --outHeader ' out_header ' --steps 3000 --threads 1']);
        [s1.vertices,s1.faces] = readPLY([out_header '.3000.ply']);
        writeOBJ(sph_map_filename, s1.vertices, s1.faces);
        delete([tempheader '/*']);
        [status, message, messageid] = rmdir(tempheader);
        if status == 0
            disp(message)
            disp(messageid)
        end
    catch
        continue
    end
    toc
end