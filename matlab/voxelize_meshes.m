%% Map area distortion to peaks
a = nd;
b = nd;
a(:)=0;
b(nd<0)=0;
for pki = 1:num_peaks
    if numel(peaks(pki).peak_vertices) > 4
        a(peaks(pki).peak_vertices) = peaks(pki).type+1;
        if peaks(pki).type == 2
            b(peaks(pki).peak_vertices) = -nd(peaks(pki).peak_vertices);
        else
            b(peaks(pki).peak_vertices) = nd(peaks(pki).peak_vertices);
        end
    end
end

%% load segmentation as ref
[o_seg_img_3d, meta] = fn_nrrdread([experiment_path '/' pid '/' pid '_CT_' nid '-seg' iso '-label.nrrd']);
B=bwperim(o_seg_img_3d,6);
S = regionprops(o_seg_img_3d,'PixelList','PixelIdxList');
if size(S.PixelList,2) == 2
    S.PixelList(:,3) = 1;
end
segmented_voxels = S.PixelList.*meta.pixelspacing'+meta.spaceorigin;
segmented_voxels(:,1:2) = -segmented_voxels(:,1:2);

%% voxelize area distortion
D1 = pdist2(s.vertices,segmented_voxels);

sigma = esr/2/3;
d = normpdf(0,0,sigma);
w1 = normpdf(D1,0,sigma)/d;
tw1=sum(w1)+1;
w = (w1./tw1)';

test_3d=double(o_seg_img_3d);
test_3d(:)=0;  

test_3d(S.PixelIdxList)=w*b;
fn_nrrdwrite([experiment_path '/' pid '/' pid '_CT_' nid '-seg-ard' iso '.nrrd'], test_3d, meta);
test_3d=test_3d.*B;
fn_nrrdwrite([experiment_path '/' pid '/' pid '_CT_' nid '-seg-ard-surface' iso '.nrrd'], test_3d, meta);


%% voxelize peak classifications
test_3d(:)=0;    
for pki = 1:num_peaks
    peaks_A = peaks(pki).l_center(1:end-1,:);
    peaks_B = peaks(pki).l_center(2:end,:);
    X = segmented_voxels;
    normals = peaks_A-peaks_B;
    d = sum(normals,2);
    normals = normals(d~=0,:); peaks_A = peaks_A(d~=0,:); peaks_B = peaks_B(d~=0,:);
    for ni = 1:size(normals,1)
        D2 = abs(normals(ni,:)*(X-peaks_B(ni,:))')./sqrt(sum(normals(ni,:).^2));
        voxels = D2<max(meta.pixelspacing)*2;
        svoxels = S.PixelIdxList(voxels);
        D3=dist([s.vertices(peaks(pki).peak_vertices,:);peaks(    pki).l_center],segmented_voxels(voxels,:)');
        minD3 = min(D3);
        test_3d(svoxels(minD3<peaks(pki).width*2/3))=peaks(pki).type+1;
    end
end
fn_nrrdwrite([experiment_path '/' pid '/' pid '_CT_' nid '-seg-peaks' iso '-label.nrrd'], test_3d, meta);
test_3d=test_3d.*B;
fn_nrrdwrite([experiment_path '/' pid '/' pid '_CT_' nid '-seg-peaks-surface' iso '.nrrd'], test_3d, meta);
