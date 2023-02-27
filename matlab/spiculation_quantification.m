function [spikes, volume, area, nd, l_ArD, l_AnD] = spiculation_quantification(s, s1, attachment_indices)
    % mesh dimensions
    n = size(s.vertices,1);
    n1 = size(s1.vertices,1);
    if n ~= n1
        a = s.faces(s.faces - s1.faces > 0);
        b = s1.faces(s.faces - s1.faces > 0);
        c = unique([a,b],'rows');
        d = c(:,1)-c(:,2);
        idx = [1; find(diff(d)<0)];
        s.vertices(c(idx,1),:) = [];
        s.faces = s1.faces;
        n = size(s.vertices,1);
    end
    
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


    %% calculate volume and surface area
    [volume,area] = stlVolume(s.vertices',s.faces');
    disp(['Orignial Volume:' num2str(volume) ', Surface area:', num2str(area)]);

    if 1 % equivilant volume sphere
        [volume1,~] = stlVolume(s1.vertices',s1.faces');
        s1.vertices = v1*(volume/volume1)^(1/3)+repmat(ce1,n,1); % make it same as original volume
    end
    [volume1,area1] = stlVolume(s1.vertices',s1.faces');
    disp(['Mapped Volume:' num2str(volume1) ', Surface area:', num2str(area1)]);


    %% compute area and angle distortions
    log_area_param = area_distortions(s, s1);
    log_angle_param = angle_distortions(s, s1);

    l_ArD = real(log_area_param(:));
    l_ArD(isinf(l_ArD)&l_ArD<0) = min(l_ArD(~isinf(l_ArD)));
    l_ArD(isinf(l_ArD)&l_ArD>0) = max(l_ArD(~isinf(l_ArD)));
    l_AnD = real(log_angle_param(:));
    l_AnD(isinf(l_AnD)&l_AnD<0) = min(l_AnD(~isinf(l_AnD)));
    l_AnD(isinf(l_AnD)&l_AnD>0) = max(l_AnD(~isinf(l_AnD)));


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


    %% spike detection
    [spikes] = detect_baselines(s, nd);
    [spikes] = detect_spikes(s, s1, Ne, nd, spikes);
    

    %% measuring heights and widths of the spikes
    [spikes] = measure_spikes(spikes,s,nd,l_ArD,attachment_indices);
end
