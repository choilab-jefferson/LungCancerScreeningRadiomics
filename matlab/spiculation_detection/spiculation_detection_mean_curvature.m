function spikes = spiculation_detection_mean_curvature(Cmean, Cgaussian, s, s1, nd, usethird)
    if ~exist('usethird','var')
        usethird = false;
    end
    
    Ne = vertex_neighbours(s); % generate neighbourhood map
    spike_candidates = find(Cmean>0.5 & abs(Cgaussian)>0.5);
    spikes = spike_candidates;
    medCmean = median(Cmean);
    while(true)
        n_spikes = numel(spikes);
        for pki = 1:numel(spike_candidates) % check negibourhood
            ppi = spike_candidates(pki);
            if(~usethird)
               Nce=unique(cell2mat(Ne(Ne{ppi})));
            else
               % Get first, second and third ring neighbours
               Nce=unique(cell2mat(Ne(cell2mat(Ne(Ne{ppi})))));
            end
            [mc, mi] = max(Cmean(Nce));
            if medCmean > 0.5 || medCmean - mc > -0.5
                spikes(pki) = 0;
                display([ppi,mc,medCmean,medCmean-mc])
            else
                spikes(pki) = Nce(mi);
            end
        end
        spikes = unique(spikes);
        spikes(spikes == 0) = [];
        n_spikes1 = numel(spikes);
        spike_candidates = spikes;
        if n_spikes == n_spikes1 || numel(spikes) == 0
            break
        end
    end
    n_spikes = numel(spikes);
    spikes_new = [];
    for pki = 1:n_spikes % check negibourhood
        spikes_new(pki).apex = spikes(pki);
        spikes_new(pki).width = 0;
        spikes_new(pki).height = 0;
        spikes_new(pki).angle = 0;
    end
    spikes = spikes_new;
	if n_spikes > 0
        sidx = [spikes.apex] > 0;
        spikes_old = spikes(sidx);
    end
    
    %% !!! area distortion baseline detecction!!!
    [spikes, t] = baseline_detection(spikes, s1, nd, Ne); 
    
    for pki = 1:n_spikes % check negibourhood
        if numel(spikes(pki).baseline) < 3
            spikes(pki).apex
        else
            [centerLoc, circleNormal, ellipse_t] = EllipseFit3D(s.vertices(spikes(pki).baseline,:));
        
            if numel(ellipse_t) == 0 || numel(ellipse_t.a) == 0
                spikes(pki).apex = 0;
            else
                [rotated_ellipse] = ellipse3D(centerLoc, circleNormal, ellipse_t);
                spikes(pki).baseline_curve = rotated_ellipse';
                spikes(pki).line_spike_base = [s.vertices(spikes(pki).apex,:);centerLoc];
                spikes(pki).height = pdist2(s.vertices(spikes(pki).apex,:),centerLoc);
                spikes(pki).width = max(pdist(s.vertices(spikes(pki).baseline,:)));

                %% angle of the spike
                nb = numel(spikes(pki).baseline);
                st.vertices = [s.vertices(spikes(pki).baseline,:);mean(s.vertices(spikes(pki).baseline,:),1);s.vertices(spikes(pki).apex,:)];
                st.faces = [(1:nb)',repmat(nb+1,nb,1),repmat(nb+2,nb,1)];
                fanglet=zeros(nb,3);
                for fi = 1:nb
                    fanglet(fi,:) = triangle_angles(st.vertices(st.faces(fi,:),:));
                end
                fanglet(isnan(fanglet)) = [];
                try
                    spikes(pki).angle=mean(fanglet(:,3)*2);
                catch
                    spikes(pki).angle=0;
                end
            end
        end
    end
    
    if n_spikes > 0
        sidx = [spikes.apex] > 0;
        spikes = spikes(sidx);
    end