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


%% voxelize spike classifications
test_3d(:)=0;    
for pki = 1:num_spikes
    spikes_A = spikes(pki).l_center(1:end-1,:);
    spikes_B = spikes(pki).l_center(2:end,:);
    X = segmented_voxels;
    normals = spikes_A-spikes_B;
    d = sum(normals,2);
    normals = normals(d~=0,:); spikes_A = spikes_A(d~=0,:); spikes_B = spikes_B(d~=0,:);
    for ni = 1:size(normals,1)
        D2 = abs(normals(ni,:)*(X-spikes_B(ni,:))')./sqrt(sum(normals(ni,:).^2));
        voxels = D2<max(meta.pixelspacing)*2;
        svoxels = S.PixelIdxList(voxels);
        D3=dist([s.vertices(spikes(pki).spike_vertices,:);spikes(    pki).l_center],segmented_voxels(voxels,:)');
        minD3 = min(D3);
        test_3d(svoxels(minD3<spikes(pki).width*2/3))=spikes(pki).type+1;
    end
end
fn_nrrdwrite([experiment_path '/' pid '/' pid '_CT_' nid '-seg-spikes' iso '-label.nrrd'], test_3d, meta);
test_3d=test_3d.*B;
fn_nrrdwrite([experiment_path '/' pid '/' pid '_CT_' nid '-seg-spikes-surface' iso '.nrrd'], test_3d, meta);
