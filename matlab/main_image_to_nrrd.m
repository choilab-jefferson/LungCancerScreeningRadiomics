clc;
clear;

%% LIDC-IDRI dicom to nrrd for radiomics
%
%  This code prepares initial data for lung cancer screening radiomics
%

set_environment; % call environment setting

% set values
%selected = selected_72;
experiment_set  = 'nodule-lidc';
iso_px_size=1; % a standard unit ('mm-unit')

% saved data load or not
load_input = true;

% directory paths
ct_img_path=[data_path '/' experiment_set '/CT'];
nodule_img_path=[data_path '/' experiment_set];

% make directory
if ~isdir(output_path); mkdir(output_path); end
if ~isdir(ct_img_path); mkdir(ct_img_path); end


%% get pids
tic
filename_pid_list = [output_path '/dicom_pid_list.mat'];
if(fn_check_load_data(filename_pid_list, load_input))
    [dicom_path_list,pid_list]=fn_scan_pid(dicom_path);

    save(filename_pid_list, 'dicom_path_list', 'pid_list');
else
    load(filename_pid_list);
end
fprintf('dicom list loaded ... \t\t\t %6.2f sec\n', toc);

all_parameters = table();

%% main process
for idx = 1:numel(pid_list)
    pid = pid_list{idx};

    %% individual selection
    %if numel(strfind(pid, '0315'))==0, continue; end

    %% selection for 72 cases
    %if sum(strcmp(selected, pid)) == 0
    %    continue
    %end

    tic % tic starts a stopwatch timer
    fprintf('%d %s\n', idx, pid);
    parameters = image_to_nrrd(ct_img_path, nodule_img_path, dicom_path_list{idx}, pid);

    all_parameters = [all_parameters; parameters];

end


% module rmpath
% rmpath('./io');
% rmpath('./interpolation');
% rmpath('./segmentation');
% rmpath('./nodule_candidate_detection');
