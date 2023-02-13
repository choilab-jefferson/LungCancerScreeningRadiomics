function s = generate_mesh_model(model_file_path, label_file_path, smooth)
    try
        [o_nodule_img_3d, meta_nodule_img_3d] = fn_nrrdread(label_file_path);
    catch exception
        throw(exception)
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

    s.vertices = (s.vertices-1).*repmat(meta.pixelspacing',size(s.vertices,1),1) + repmat(meta.spaceorigin,size(s.vertices,1),1);
    if meta.spacedirections(3,3) < 0
        s.vertices(:,1:2) = -s.vertices(:,1:2);
    end
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

    s.vertices = v;
    s.faces = f;

    writeOBJ(model_file_path, v, f);
end