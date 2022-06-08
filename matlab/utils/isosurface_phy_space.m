function [s] = isosurface_phy_space(BW, meta)
    smooth = true;
    s = isosurface(BW, 0.4);
    n = size(s.vertices,1);
    m = size(s.faces,1);
    s.vertices = s.vertices.*repmat(meta.pixelspacing',n,1) + repmat(meta.spaceorigin,n,1);
    s.vertices(:,1:2) = -s.vertices(:,1:2);
    if smooth
        s = smoothpatch(s,1,1);
    end
    v = s.vertices; f = s.faces;
    try
        [v,f]=meshcheckrepair(v,f); % iso2mesh
    catch
        [v,f]=meshcheckrepair(v,f,'meshfix');
    end
    %[v,f]=remeshsurf(v,f,1);
    [v,f]=surfreorient(v,f);
    [v,f]=meshcheckrepair(v,f);
end
