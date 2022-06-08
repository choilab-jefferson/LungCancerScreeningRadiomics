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
        rnid = num2str(nodule_info{nid, 'nid'});
        try
            [o_nodule_img_3d, meta_nodule_img_3d] = fn_nrrdread([experiment_path '/' pid '/' pid '_CT_' rnid '-seg' iso '-label.nrrd']);
        catch
            continue
        end
        sz_o_nodule_img_3d = size(o_nodule_img_3d);

        meta = struct();
        meta.type = 'float';
        meta.encoding = 'gzip';
        meta.endian = 'little';
        meta.spacedirections = meta_nodule_img_3d.spacedirections;
        meta.pixelspacing = meta_nodule_img_3d.pixelspacing;
        
        
        sz_nodule_blk_3d = max([sz_bbox,sz_bbox,sz_bbox],sz_o_nodule_img_3d);
        offsetXYZ = round((sz_nodule_blk_3d - sz_o_nodule_img_3d)/2)+1;
        BW = zeros(sz_nodule_blk_3d);
        BW(offsetXYZ(2):offsetXYZ(2)+sz_o_nodule_img_3d(1)-1, ...
            offsetXYZ(1):offsetXYZ(1)+sz_o_nodule_img_3d(2)-1, ...
            offsetXYZ(3):offsetXYZ(3)+sz_o_nodule_img_3d(3)-1) = o_nodule_img_3d;
        display(size(BW))
        
        meta.spaceorigin = meta_nodule_img_3d.spaceorigin-offsetXYZ.*meta.pixelspacing';
        
        
        %% select the biggest object
        ce_bbox = round(sz_bbox/2);
        CC = bwconncomp(BW,6);
        idx1 = find(cellfun(@(x) sum(x == sz_bbox*sz_bbox*(ce_bbox-1)+sz_bbox*(ce_bbox-1)+ce_bbox), CC.PixelIdxList));
        if numel(idx1) == 0
            numPixels = cellfun(@numel,CC.PixelIdxList);
            [biggest, idx1] = max(numPixels);
        end
        BW(:) = 0;
        BW(CC.PixelIdxList{idx1}) = 1;

        %%
        if sum(BW(:)) == 0
            continue
        end


        merged_nodule_info = [merged_nodule_info;nodule_info(nid,:)];
        
        %%
        S = regionprops(BW,'Centroid');
        ce = S.Centroid;
        display(ce)
             
        if ~smooth
            s = isosurface(BW,0.4);
        else
            %Lempisky smooth surface optimization
            V = logical(BW);
            [f,v] = ExtractSeparatingSurfaceCPU(V, 1000);
            s.vertices = v; s.faces = f;
        end
        n = size(s.vertices,1);
        m = size(s.faces,1);
        s.vertices = s.vertices.*repmat(meta.pixelspacing',n,1) + repmat(meta.spaceorigin,n,1);
        s.vertices(:,1:2) = -s.vertices(:,1:2);

        v = s.vertices; f = s.faces;
        try
            [v,f]=meshcheckrepair(v,f); % iso2mesh
        catch
            [v,f]=meshcheckrepair(v,f,'meshfix');
        end
        %[v,f]=remeshsurf(v,f,1);
        [v,f]=surfreorient(v,f);
        [v,f]=meshcheckrepair(v,f);
        try
            statistics(v,f) % gptoolbox
        catch
        end
        writeOBJ([objs_experiment_path '/' pid '_' rnid '.obj'], v, f);
        %break
    end
    close all

    fclose all;
    toc
end

writetable(merged_nodule_info, [objs_experiment_path '/nodule_info.csv'])

