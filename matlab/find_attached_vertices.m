function attachment = find_attached_vertices(o_attached_img_3d, meta, s)
    [x,y,z] = ndgrid(-1:1);
    se1 = strel(sqrt(x.^2 + y.^2 + z.^2) <=1);

    try
        attached_img_3d = imclose(imerode(o_attached_img_3d,se1),se1);
        S = regionprops(attached_img_3d,'PixelList');
        nn = size(S.PixelList,1);
        if size(S,1) > 0 && nn>1
            attached_voxels = S.PixelList.*repmat(meta.pixelspacing',nn,1)+repmat(meta.spaceorigin,nn,1);
            if meta.spacedirections(3,3) < 0
                attached_voxels(:,1:2) = -attached_voxels(:,1:2);
            end

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