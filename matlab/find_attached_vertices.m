function attachment = find_attached_vertices(o_attached_img_3d, meta, s)
    [x,y,z] = ndgrid(-1:1);
    se1 = strel(sqrt(x.^2 + y.^2 + z.^2) <=1);

    try
        attached_img_3d = imclose(imerode(o_attached_img_3d,se1),se1);
        S = regionprops(bwperim(attached_img_3d, 26),'PixelList');
        pixelList = vertcat(S.PixelList);
        nn = size(pixelList,1);
        if size(pixelList,1) > 0 && nn>1
            attached_voxels = (pixelList-1).*repmat(diag(meta.spacedirections)',nn,1) + repmat(meta.spaceorigin,nn,1);

            D = pdist2(attached_voxels,s.vertices);
            attachment = (min(D)<meta.pixelspacing(1)*2)';
        else
            attachment=zeros(size(s.vertices,1),1);
        end
    catch exception
        disp(exception)
        attachment=zeros(size(s.vertices,1),1);
    end
end