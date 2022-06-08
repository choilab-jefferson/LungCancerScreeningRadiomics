if ~exist('refine','var')
    refine = false;
end
if ~exist('smooth','var')
    smooth = false;
end

th_noise = 3;

results = table();
tic % tic starts a stopwatch timer
try
    nodule_info = readtable([experiment_path '/' pid '/' pid '_all.csv']);
    if numel(nids) == 1
        n_info = nodule_info(nodule_info.nid == str2double(nid),:);
    else % use original nid if there are multiple segmentations
        n_info = nodule_info(nodule_info.nid == str2double(nids{1}),:);
    end
catch
    n_info = [];        
end


%% load nodule mesh model and its spherical mapping
try
    if ~isempty(nid)
        obj_filename = [obj_path '/' pid '_' nid '.obj'];
        sph_map_filename = [sph_map_path '/' pid '_' nid '_spherical.obj'];
    else
        obj_filename = [obj_path '/' pid '.obj'];
        sph_map_filename = [sph_map_path '/' pid '_spherical.obj'];
    end

    [s.vertices,s.faces] = readOBJ(obj_filename);
    if refine, s = refinepatch(s); end
    if smooth, s = smoothpatch(s,1,1); end
    try
        [s1.vertices,s1.faces,~,normals] = read_obj_colors(sph_map_filename);
    catch
        [s1.vertices,s1.faces] = readOBJ(sph_map_filename);
        [normals,normals_f] = patchnormals(s1);
    end
    if refine, s1 = refinepatch(s1); end
    if smooth, s1 = smoothpatch(s1,1,1); end
    assert(numel(s.vertices)==numel(s1.vertices));
catch exception
    if ~exist(obj_filename, 'file')
        display('no obj file')
    end
    if ~exist(sph_map_filename, 'file') || numel(s.vertices)~=numel(s1.vertices)
        display('no spherical mapping!')
        s1=[];
        normals=[];
        %[s1,~,normals] = sphereical_mapping(s);
    end
    
end

%%

% mesh dimensions
n = size(s.vertices,1);
m = size(s.faces,1);
ce = mean(s.vertices);

v = s.vertices-repmat(ce,n,1); % centering the nodule mesh model
%s.vertices = v;
[normal_s,normalf_s] = patchnormals(s);
colors = (1-normal_s)/2; % recaculate vetex color
colors = 1-colors; % flip color space to show bright color

% normalization of the model
magv = sqrt(sum(v.^2,2));
mdist = mean(magv);
nv = v./repmat(magv,1,3); % nodule

ce1 = mean(s1.vertices);
v1 = s1.vertices-repmat(ce1,n,1);
magv1 = sqrt(sum(v1.^2,2));
nv1 = v1./repmat(magv1,1,3); % spherical mapping
%s1.vertices = nv1*mdist+repmat(ce1,n,1);

Ne=vertex_neighbours(s); % generate neighbourhood map


%% load attachment information
try
    [o_attached_img_3d, meta] = fn_nrrdread([experiment_path '/' pid '/' pid '_CT_' nid '-attached-label.nrrd']);
    attached_img_3d = imclose(imerode(o_attached_img_3d,se1),se1);
    S = regionprops(attached_img_3d,'PixelList');
    nn = size(S.PixelList,1);
    if size(S,1) > 0 && nn>1
        attached_voxels = S.PixelList.*repmat(meta.pixelspacing',nn,1)+repmat(meta.spaceorigin,nn,1);
        attached_voxels(:,1:2) = -attached_voxels(:,1:2);
        
        D = pdist2(attached_voxels,s.vertices);
        attachment = (min(D)<meta.pixelspacing(1)*2)';
        attached_vertices=find(attachment);
    else
        attachment=zeros(size(s.vertices,1),1);
        attached_vertices=[];
    end
catch
    attachment=zeros(size(s.vertices,1),1);
    attached_vertices=[];
end


%% calculate volume and surface area
[volume,area] = stlVolume(s.vertices',s.faces');
disp(['Orignial Volume:' num2str(volume) ', Surface area:', num2str(area)]);

if 1 % equivilant volume sphere
    [volume1,~] = stlVolume(s1.vertices',s1.faces');
    s1.vertices = v1*(volume/volume1)^(1/3)+repmat(ce1,n,1); % make it same as original volume
end
[volume1,area1] = stlVolume(s1.vertices',s1.faces');
disp(['Mapped Volume:' num2str(volume1) ', Surface area:', num2str(area1)]);


%% extract shape features from nodule mesh model
esr = esph_factor*volume^(1/3);
roundness = min((4*pi*esr^2)/area,1);
roundness1 = min(area1/area,1);


%% compute area and angle distortions
log_area_param = area_distortions(s, s1);
log_angle_param = angle_distortions(s, s1);

l_ArD = real(log_area_param(:));l_ArD(isinf(l_ArD)&l_ArD<0) = min(l_ArD(~isinf(l_ArD)));l_ArD(isinf(l_ArD)&l_ArD>0) = max(l_ArD(~isinf(l_ArD)));
l_AnD = real(log_angle_param(:));l_AnD(isinf(l_AnD)&l_AnD<0) = min(l_AnD(~isinf(l_AnD)));l_AnD(isinf(l_AnD)&l_AnD>0) = max(l_AnD(~isinf(l_AnD)));


%% raw features from curvature maps
[Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2]=patchcurvature(s);

l_GC = Cgaussian(:);
l_MC = Cmean(:); 


%%
[muhat,sigmahat] = normfit(l_ArD); % Gaussian fitting
if 0 % Gaussian fitting with censoring
    pmuhat = 0; psigmahat = 0;
    censoring = zeros(size(l_ArD));
    while (abs(muhat-pmuhat) > eps || abs(sigmahat-psigmahat) > eps) && sum(censoring)/numel(censoring) < 0.1
        pmuhat = muhat; psigmahat = sigmahat;
        censoring = l_ArD < muhat-sigmahat*3; % cover 99.87%
        [muhat,sigmahat] = normfit(l_ArD,[],censoring); 
    end
end
nd = (l_ArD - muhat)/sigmahat;

%% peak detection by mean curvature
if 0
    peaks_old = spiculation_detection_mean_curvature(Cmean, Cgaussian, Ne, s, s1, nd);
end

%% peak detection
[peaks] = detect_baselines(s, nd);
[peaks] = detect_peaks(s, s1, Ne, nd, peaks);

%% measuring heights and widths of the peaks
[peaks] = measure_peaks(peaks,s,nd,l_ArD,attached_vertices);


%%
%figure, hist([peaks.height])
[num_peaks,l_height,l_width,l_angle,l_nheight,l_nwidth,l_rimchange, ...
      spic_ap,spic_bp,spic1p,spic2p,spic3p,spic4p,hx,hn] = extract_spiculation_features(peaks, esr);
nhn = hn / sum(hn);

%% convert peaks structure to table
if num_peaks > 0
    if num_peaks > 1
        peaks_table = struct2table(peaks);
    else
        peaks(2) = peaks(1);
        peaks_table = struct2table(peaks);
        peaks_table = peaks_table(1,:);
        peaks(2) = [];
    end
    peaks_table.PID = repmat({pid},num_peaks,1);
    peaks_table.NID = repmat({nid},num_peaks,1);
    peaks_table.D = repmat(esr*2,num_peaks,1);
    peaks_table.volume = repmat(volume/1000,num_peaks,1);
    peaks_table.nd = nd(peaks_table.apex);
    peaks_table.rimchange = l_rimchange';
    peaks_table.nrimchange = l_rimchange'./(l_height'+1);
    if numel(n_info) > 0
        peaks_table.spic = repmat(n_info.Spiculation,num_peaks,1);
        peaks_table.lobul = repmat(n_info.Lobulation,num_peaks,1);
        peaks_table.malig = repmat(n_info.Malignancy,num_peaks,1);
        peaks_table = peaks_table(:,{'PID','NID','spic','lobul','malig','volume','D','nd','nrimchange','apex','height','width','angle'});
    else
        peaks_table = peaks_table(:,{'PID','NID','volume','D','nd','nrimchange','apex','height','width','angle'});
    end
    
    disp(peaks_table)
    all_peaks = [all_peaks; peaks_table];
end


%% remove small peaks
selected_peaks = l_height > th_noise & l_width > th_noise/2;
peaks = peaks(selected_peaks);
num_peaks = numel(peaks);
fprintf('number of peaks %d\n', num_peaks)

%%
l_loc = s.vertices([peaks.apex],:); % list of apex coordinates


%% peak classification, spiculation or lobulation
spiculation = [peaks.type]==0;
lobulation = [peaks.type]==1;
attached = [peaks.type]==3;


%% recalculate features
[num_peaks,l_height,l_width,l_angle,l_nheight,l_nwidth,l_rimchange, ...
      spic_a,spic_b,spic1,spic2,spic3,spic4,hx1,hn1] = extract_spiculation_features(peaks, esr);
nhn1 = hn1 / sum(hn1);


%% save important parameters
dlmwrite([output_experiment_path '/parameters/' pid '_' nid '_nd.txt'],nd);
dlmwrite([output_experiment_path '/parameters/' pid '_' nid '_attached.txt'], attached_vertices);
dlmwrite([output_experiment_path '/parameters/' pid '_' nid '_apex.txt'],[peaks.apex]');

%% save apex point of peaks and area distortion for slicer
fid = fopen([output_experiment_path '/parameters/' pid '_' nid '.fcsv'],'w');
fprintf(fid, '# Markups fiducial file version = 4.8\n');
fprintf(fid, '# CoordinateSystem = 0\n');
fprintf(fid, '# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n');
for pki = 1:num_peaks
    if peaks(pki).type == 0
        fprintf(fid, 'vtkMRMLMarkupsFiducialNode_%d,%f,%f,%f,0,0,0,0,1,1,1,S%d,"AD %0.2f",vtkMRMLScalarVolumeNode1\n',pki,l_loc(pki,1),l_loc(pki,2),l_loc(pki,3),pki,nd(peaks(pki).apex));
    elseif peaks(pki).type == 1
        fprintf(fid, 'vtkMRMLMarkupsFiducialNode_%d,%f,%f,%f,0,0,0,0,1,0,2,L%d,"AD %0.2f",vtkMRMLScalarVolumeNode1\n',pki,l_loc(pki,1),l_loc(pki,2),l_loc(pki,3),pki,nd(peaks(pki).apex));
    elseif peaks(pki).type == 2
        fprintf(fid, 'vtkMRMLMarkupsFiducialNode_%d,%f,%f,%f,0,0,0,0,0,0,2,C%d,"AD %0.2f",vtkMRMLScalarVolumeNode1\n',pki,l_loc(pki,1),l_loc(pki,2),l_loc(pki,3),pki,nd(peaks(pki).apex));
    elseif peaks(pki).type == 3
        fprintf(fid, 'vtkMRMLMarkupsFiducialNode_%d,%f,%f,%f,0,0,0,0,0,0,2,A%d,"AD %0.2f",vtkMRMLScalarVolumeNode1\n',pki,l_loc(pki,1),l_loc(pki,2),l_loc(pki,3),pki,nd(peaks(pki).apex));
    end
end
fclose(fid);

%writeOBJ([output_experiment_path '/' pid '_' nid '.obj'],s.vertices,s.faces,nd,s.faces,normal_s,normalf_s);



%% feature table
hnames = strsplit(sprintf('h%d ',1:10));
nhnames = strsplit(sprintf('nh%d ',1:10));

f = table(idx,{pid},{nid},num_peaks,sum(spiculation),sum(lobulation), sum(attached), mean(attachment), ...
        volume,area,esr*2,roundness,spic_a, spic_b, spic_ap, spic_bp, ...
        spic1,spic2,spic3,spic4,spic1p,spic2p,spic3p,spic4p,...
        min(l_GC),max(l_GC),median(l_GC),mean(l_GC),var(l_GC),skewness(l_GC),kurtosis(l_GC), ...
        min(l_MC), max(l_MC),median(l_MC),mean(l_MC),var(l_MC),skewness(l_MC),kurtosis(l_MC), ...
        min(nd),max(nd),median(nd),mean(nd),var(nd),skewness(nd),kurtosis(nd), ...
        min(l_ArD),max(l_ArD),median(l_ArD),mean(l_ArD),var(l_ArD),skewness(l_ArD),kurtosis(l_ArD), ...
        min(l_AnD), max(l_AnD),median(l_AnD),mean(l_AnD),var(l_AnD),skewness(l_AnD),kurtosis(l_AnD), ...
        min(l_height), max(l_height),median(l_height),mean(l_height),var(l_height),skewness(l_height),kurtosis(l_height), ...
        min(l_width), max(l_width),median(l_width),mean(l_width),var(l_width),skewness(l_width),kurtosis(l_width), ...
        min(l_angle), max(l_angle),median(l_angle),mean(l_angle),var(l_angle),skewness(l_angle),kurtosis(l_angle), ...
        min(l_nheight), max(l_nheight),median(l_nheight),mean(l_nheight),var(l_nheight),skewness(l_nheight),kurtosis(l_nheight), ...
        min(l_nwidth), max(l_nwidth),median(l_nwidth),mean(l_nwidth),var(l_nwidth),skewness(l_nwidth),kurtosis(l_nwidth), ...
        hn(1),hn(2),hn(3),hn(4),hn(5),hn(6),hn(7),hn(8),hn(9),hn(10), ...
        nhn(1),nhn(2),nhn(3),nhn(4),nhn(5),nhn(6),nhn(7),nhn(8),nhn(9),nhn(10), ...
        'VariableNames',{'No','PID','NID','num_peaks','num_spic','num_lobul','num_attached','attachment', ...
        'volume','area','eqv_sph_D','roundness','spic_a','spic_b','spic_ap','spic_bp', ...
        'spic1','spic2','spic3','spic4','spic1p','spic2p','spic3p','spic4p', ...
        'min_GC','max_GC','med_GC','mean_GC','var_GC','skew_GC','kurt_GC', ...
        'min_MC','max_MC','med_MC','mean_MC','var_MC','skew_MC','kurt_MC', ...
        'min_NArD','max_NArD','med_NArD','mean_NArD','var_NArD','skew_NArD','kurt_NArD', ...
        'min_ArD','max_ArD','med_ArD','mean_ArD','var_ArD','skew_ArD','kurt_ArD', ...
        'min_AnD','max_AnD','med_AnD','mean_AnD','var_AnD','skew_AnD','kurt_AnD', ...
        'min_h','max_h','med_h','mean_h','var_h','skew_h','kurt_h', ...
        'min_w','max_w','med_w','mean_w','var_w','skew_w','kurt_w', ...
        'min_a','max_a','med_a','mean_a','var_a','skew_a','kurt_a', ...
        'min_nh','max_nh','med_nh','mean_nh','var_nh','skew_nh','kurt_nh', ...
        'min_nw','max_nw','med_nw','mean_nw','var_nw','skew_nw','kurt_nw' ...
         hnames{1:10},nhnames{1:10}});


display(f)
features = [features;f];