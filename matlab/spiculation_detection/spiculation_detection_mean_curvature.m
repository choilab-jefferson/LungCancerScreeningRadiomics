function peaks = spiculation_detection_mean_curvature(Cmean, Cgaussian, Ne, s, s1, nd, usethird)
    if ~exist('usethird','var')
        usethird = false;
    end
    
    peak_candidates = find(Cmean>0.5 & abs(Cgaussian)>0.5);
    peaks = peak_candidates;
    medCmean = median(Cmean);
    while(true)
        n_peaks = numel(peaks);
        for pki = 1:numel(peak_candidates) % check negibourhood
            ppi = peak_candidates(pki);
            if(~usethird)
               Nce=unique(cell2mat(Ne(Ne{ppi})));
            else
               % Get first, second and third ring neighbours
               Nce=unique(cell2mat(Ne(cell2mat(Ne(Ne{ppi})))));
            end
            [mc, mi] = max(Cmean(Nce));
            if medCmean > 0.5 || medCmean - mc > -0.5
                peaks(pki) = 0;
                display([ppi,mc,medCmean,medCmean-mc])
            else
                peaks(pki) = Nce(mi);
            end
        end
        peaks = unique(peaks);
        peaks(peaks == 0) = [];
        n_peaks1 = numel(peaks);
        peak_candidates = peaks;
        if n_peaks == n_peaks1 || numel(peaks) == 0
            break
        end
    end
    n_peaks = numel(peaks);
    peaks_new = [];
    for pki = 1:n_peaks % check negibourhood
        peaks_new(pki).apex = peaks(pki);
        peaks_new(pki).width = 0;
        peaks_new(pki).height = 0;
        peaks_new(pki).angle = 0;
    end
    peaks = peaks_new;
	if n_peaks > 0
        sidx = [peaks.apex] > 0;
        peaks_old = peaks(sidx);
    end
    
    %% !!! area distortion baseline detecction!!!
    [peaks, t] = baseline_detection(peaks, s1, nd, Ne); 
    
    for pki = 1:n_peaks % check negibourhood
        if numel(peaks(pki).baseline) < 3
            peaks(pki).apex
        else
            [centerLoc, circleNormal, ellipse_t] = EllipseFit3D(s.vertices(peaks(pki).baseline,:));
        
            if numel(ellipse_t) == 0 || numel(ellipse_t.a) == 0
                peaks(pki).apex = 0;
            else
                [rotated_ellipse] = ellipse3D(centerLoc, circleNormal, ellipse_t);
                peaks(pki).baseline_curve = rotated_ellipse';
                peaks(pki).line_peak_base = [s.vertices(peaks(pki).apex,:);centerLoc];
                peaks(pki).height = pdist2(s.vertices(peaks(pki).apex,:),centerLoc);
                peaks(pki).width = max(pdist(s.vertices(peaks(pki).baseline,:)));

                %% angle of the peak
                nb = numel(peaks(pki).baseline);
                st.vertices = [s.vertices(peaks(pki).baseline,:);mean(s.vertices(peaks(pki).baseline,:),1);s.vertices(peaks(pki).apex,:)];
                st.faces = [(1:nb)',repmat(nb+1,nb,1),repmat(nb+2,nb,1)];
                fanglet=zeros(nb,3);
                for fi = 1:nb
                    fanglet(fi,:) = triangle_angles(st.vertices(st.faces(fi,:),:));
                end
                fanglet(isnan(fanglet)) = [];
                try
                    peaks(pki).angle=mean(fanglet(:,3)*2);
                catch
                    peaks(pki).angle=0;
                end
            end
        end
    end
    
    if n_peaks > 0
        sidx = [peaks.apex] > 0;
        peaks = peaks(sidx);
    end