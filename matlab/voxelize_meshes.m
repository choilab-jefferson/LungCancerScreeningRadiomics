%% Map area distortion to spikes
a = nd;
b = nd;
a(:)=0;
b(nd<0)=0;
for pki = 1:num_spikes
    if numel(spikes(pki).spike_vertices) > 4
        a(spikes(pki).spike_vertices) = spikes(pki).type+1;
        if spikes(pki).type == 2
            b(spikes(pki).spike_vertices) = -nd(spikes(pki).spike_vertices);
        else
            b(spikes(pki).spike_vertices) = nd(spikes(pki).spike_vertices);
        end
    end
end

%% load segmentation as ref
label_file = dir([experiment_path '/' pid '/' pid '_CT_' nid '-R1' iso_t '-label.nrrd']);
[o_seg_img_3d, meta] = fn_nrrdread([label_file.folder '/' label_file.name]);
B=bwperim(o_seg_img_3d,6);
S = regionprops(o_seg_img_3d,'PixelList','PixelIdxList');
if size(S.PixelList,2) == 2
    S.PixelList(:,3) = 1;
end
segmented_voxels = S.PixelList.*meta.pixelspacing'+meta.spaceorigin;
segmented_voxels(:,1:2) = -segmented_voxels(:,1:2);

%% voxelize area distortion
sigma = esr/2/3;
d = normpdf(0,0,sigma);
test_3d=double(o_seg_img_3d);
test_3d(:)=0;  

if size(segmented_voxels, 1) * size(s.vertices, 1) < 1024*1024*1024*10
    D1 = pdist2(s.vertices,segmented_voxels);
    w1 = normpdf(D1,0,sigma)/d;
    tw1 = sum(w1);
    tw1(tw1 == 0) = 1;
    w = (w1./tw1)';

    test_3d(S.PixelIdxList)=w*b;
else
    t = zeros(1, size(segmented_voxels, 1));
    parfor vidx = 1:size(segmented_voxels, 1)
        v = segmented_voxels(vidx, :);
        D1 = pdist2(s.vertices, v);
        w1 = normpdf(D1,0,sigma)/d;
        %w1(w1 < 1e-10) = 0;
        tw1 = sum(w1);
        if tw1 == 0
            tw1 = 1
        end
        w = (w1./tw1)';
        t(vidx) = w*b;
    end
    test_3d(S.PixelIdxList)=t;
end

ard_3d = test_3d;
fn_nrrdwrite([experiment_path '/' pid '/' pid '_CT_' nid '-R1-ard' iso_t smooth_t '.nrrd'], test_3d, meta);
test_3d=test_3d.*B;
fn_nrrdwrite([experiment_path '/' pid '/' pid '_CT_' nid '-R1-ard-surface' iso_t smooth_t '.nrrd'], test_3d, meta);


%% voxelize spike classifications
spikes_voxel_idx = find(ard_3d(S.PixelIdxList) < 0);
test_3d(:)=0;    
for pki = 1:num_spikes
    fv = s;
    fv.faces = s.faces(spikes(pki).faces,:);
    fv = reducepatch(fv);
    D = pdist2(spikes(pki).l_center, segmented_voxels(spikes_voxel_idx,:));
    minD = min(D);
    selected_voxel_idx = find(minD < spikes(pki).width*2);
    X = segmented_voxels(spikes_voxel_idx(selected_voxel_idx), :);
    
    spikes_B = spikes(pki).l_center(2:end,:);
    normals = spikes(pki).l_normal;
    max_d = max(max(dist(spikes(pki).l_center')));
    for ni = 1:size(normals,1)
        D2 = abs(normals(ni,:)*(X-spikes_B(ni,:))')./sqrt(sum(normals(ni,:).^2));
        voxels = D2<max_d*3;
        svoxels = S.PixelIdxList(spikes_voxel_idx(selected_voxel_idx(voxels)));
        %D3=dist([s.vertices(spikes(pki).spike_vertices,:);spikes(pki).l_center],X(voxels,:)');
        D3=dist([fv.vertices;spikes(pki).l_center],X(voxels,:)');
        minD3 = min(D3);
        test_3d(svoxels(minD3<spikes(pki).width*2/3))=spikes(pki).type+1;
    end
end
fn_nrrdwrite([experiment_path '/' pid '/' pid '_CT_' nid '-R1-spikes' iso_t smooth_t '-label.nrrd'], test_3d, meta);
test_3d=test_3d.*B;
fn_nrrdwrite([experiment_path '/' pid '/' pid '_CT_' nid '-R1-spikes-surface' iso_t smooth_t '-smooth.nrrd'], test_3d, meta);
