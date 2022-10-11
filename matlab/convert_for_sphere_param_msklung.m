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
    label_files = dir([experiment_path '/' pid '/' pid '_CT_*-R1' iso '-label.nrrd']);
    parfor l = 1:numel(label_files)
        tic
        label_file = label_files(l);
        %% input part
        label_file_path = [label_file.folder '/' label_file.name];
        try
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
        id = strsplit(label_file.name,'_');
        pid = id{1};
        id = strsplit(id{3},'-');
        nid = [id{1} '-' id{2}];
%         sum(attached_vertices)/size(s.vertices,1)
%         patch(s, 'FaceVertexCData',double(attached_vertices),'FaceAlpha',1,'FaceColor','interp', 'EdgeColor', 'none'),caxis([0, 1]),view([-37.5 30]);
%         print(gcf,[objs_experiment_path '/' pid '_' nid '_attachement.png'],'-dpng','-r300')
        writeOBJ([objs_experiment_path '/' pid '_' nid '.obj'], v, f);
        %break

        close all

        fclose all;
        toc
    end
end

%writetable(merged_nodule_info, 'nodule_info.csv')

