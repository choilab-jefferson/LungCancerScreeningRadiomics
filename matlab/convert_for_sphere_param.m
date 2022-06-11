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
    %% input part
    nodule_info = readtable([experiment_path '/' pid '/' pid '_all.csv']);
    nodule_info(isnan(nodule_info.Volume),:) =[];
    nodule_info = sortrows(nodule_info, 'Volume', 'descend');

    if ~isdir(['output/' pid]); mkdir(['output/' pid]); end
    for nid = 1:size(nodule_info,1)
        %% input part
        rnid = num2str(nodule_info{nid, 'nid'});
        try
            label_file_path = [experiment_path '/' pid '/' pid '_CT_' rnid '-seg' iso '-label.nrrd'];
            [o_nodule_img_3d, meta_nodule_img_3d] = fn_nrrdread(label_file_path);
%             [o_attached_img_3d, ~] = fn_nrrdread(strrep(label_file_path, ['seg' iso], 'attached'));
        catch
            continue
        end
        sz_o_nodule_img_3d = size(o_nodule_img_3d);
        if numel(sz_o_nodule_img_3d) == 2
            sz_o_nodule_img_3d(3) = 1;
        end

        meta = struct();
        meta.type = 'float';
        meta.encoding = 'gzip';
        meta.endian = 'little';
        meta.spacedirections = meta_nodule_img_3d.spacedirections;
        meta.pixelspacing = meta_nodule_img_3d.pixelspacing;
        meta.spaceorigin = meta_nodule_img_3d.spaceorigin;

        BW = o_nodule_img_3d;


        %%
        S = regionprops(BW>0.4,'Centroid');
        ce = S.Centroid;
        display(ce)

        s = isosurface(BW,0.4);
        n = size(s.vertices,1);
        m = size(s.faces,1);
%         S = regionprops(imerode(o_attached_img_3d,se1),'PixelList');
%         if size(S,1) > 0 && size(S.PixelList,1)>1
%             D = pdist2(S.PixelList,s.vertices);
%             attached_vertices=(min(D)<3)';   
%         else
%             attached_vertices=zeros(size(s.vertices,1),1);
%         end
        s.vertices = (s.vertices-1).*repmat(meta.pixelspacing',size(s.vertices,1),1) + repmat(meta.spaceorigin,size(s.vertices,1),1);
        s.vertices(:,1:2) = -s.vertices(:,1:2);
        if smooth
            s = smoothpatch(s,1,1);
        end
        v = s.vertices; f = s.faces;
        [v,f]=meshcheckrepair(v,f); % iso2mesh
        [v,f]=meshcheckrepair(v,f,'meshfix');
        %[v,f]=remeshsurf(v,f,1);
        [v,f]=surfreorient(v,f);
        [v,f]=meshcheckrepair(v,f);
        try
            statistics(v,f) % gptoolbox
        catch
        end
%         sum(attached_vertices)/size(s.vertices,1)
%         patch(s, 'FaceVertexCData',double(attached_vertices),'FaceAlpha',1,'FaceColor','interp', 'EdgeColor', 'none'),caxis([0, 1]),view([-37.5 30]);
%         print(gcf,[objs_experiment_path '/' pid '_' rnid '_attachement.png'],'-dpng','-r300')
        writeOBJ([objs_experiment_path '/' pid '_' rnid '.obj'], v, f);
        %break

        close all

        fclose all;
        toc
    end
end

writetable(merged_nodule_info, [objs_experiment_path '/nodule_info.csv'])

