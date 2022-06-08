function [peaks] = detect_peaks(s, s1, Ne, nd, peaks)
%     patch(s, 'FaceVertexCData',nd,'FaceAlpha',0.4,'FaceColor','interp', 'EdgeColor', 'none');
%     colormap('jet'), caxis([-1 1])
%     view(-37.5,30);
%     axis equal fill
%     hold on

    global geodesic_library;                
    geodesic_library = 'geodesic_release';      %"release" is faster and "debug" does additional checks
    
    for i = 1:numel(peaks)
        peaks = find_next_rim(s, s1, Ne, nd, {peaks(i).baseline}, peaks, i);
    end
    idx = [peaks.apex]>0;
    peaks = peaks(idx);
    
    geodesic_delete;
end

function [peaks] = find_next_rim(s, s1, Ne, nd, rims, peaks, parent)
    nr = numel(rims); % number of rims
    
    %% calculate perimeter of rims
    p_rims = zeros(nr,1);
    for rid = 1:nr
        rim = rims{rid};
        v_rim = s.vertices(rim,:);
        p_rims(rid) = sum(dists_points(v_rim([1:end,1],:))); 
    end
    
    %% find next rim
    for rid = 1:nr
        rim = rims{rid};
        nrims = {};
        %disp([numel(rims),numel(rim)])
        %plot3(s.vertices(rim,1),s.vertices(rim,2),s.vertices(rim,3),'r')
        
        v_rim = s.vertices(rim,:);
        [center,normal] = calcaulte_center_normal(v_rim);
        
        pki = numel(peaks)+1;
        peaks(pki) = peaks(parent);
        peaks(pki).parent = parent;
        
        % when there are multiple rims and its perimeter is less than 60% of entire perimeter
        % replace baseline
        if 1 && nr > 1 && p_rims(rid)/sum(p_rims) < 0.6 && dists_points([peaks(pki).center; center]) > 2*peaks(parent).height/numel(peaks(parent).rims)
            peaks(pki).baseline = peaks(parent).rims{end};
            peaks(pki).peak_vertices = peaks(parent).vertices;
            peaks(pki).baseline_curve = s.vertices(rim([1:end,1]),:);
            peaks(pki).rims = {rim};
            peaks(pki).l_center = center;
            peaks(pki).l_normal = normal;
            fprintf('baseline %d %d %d %0.2f\n',pki, parent, nr, p_rims(rid))
        else % add a rim into the peak
            peaks(pki).rims = [peaks(pki).rims;rim];
            peaks(pki).l_center = [peaks(pki).l_center;center];
            peaks(pki).l_normal = [peaks(pki).l_normal;normal];
        end

        % update peak measures
        peaks(pki).height = sum(dists_points(peaks(pki).l_center));
        peaks(pki).center = center;
        peaks(pki).normal = normal;
                
        
        %% find inside vertices
        v_rim1 = s1.vertices(rim,:);
        v_peak1 = s1.vertices(peaks(pki).vertices,:);
        [~,normal1] = calcaulte_center_normal(v_rim1);
        xy_rim = RodriguesRotation(v_rim1,normal1,[0, 0, 1]);
        xy_peak = RodriguesRotation(v_peak1,normal1,[0, 0, 1]);
        dist2rim = min(pdist2(v_rim1,v_peak1));
        in_vertices = inpolygon(xy_peak(:,1),xy_peak(:,2),xy_rim(:,1),xy_rim(:,2));
        in_vertices1 = (dist2rim>0&dist2rim<max(pdist(v_rim1)))'; % not in the rim
                
        peaks(pki).vertices = peaks(pki).vertices(in_vertices&in_vertices1); % select inside vertices

        medndr = median(nd(rim));
        mndr = mean(nd(rim));
        %next_rim = peaks(pki).vertices(nd(peaks(pki).vertices)<min(medndr,mndr));
        next_rim = peaks(pki).vertices(nd(peaks(pki).vertices)<min(nd(rim)));

%         v_inside = s.vertices(peaks(pki).vertices,:);
%         plot3(v_inside(:,1),v_inside(:,2),v_inside(:,3),'+'), hold on,
%         plot3(v_rim([1:end,1],1),v_rim([1:end,1],2),v_rim([1:end,1],3),'*'), hold on,
%         plot3(v_rim([1:end,1],1),v_rim([1:end,1],2),v_rim([1:end,1],3)), hold on,

        %% next rim
        if numel(next_rim) > 4
            s_next = s;
            is_rim_face = sum(ismember(s.faces,next_rim),2)==3;
            s_next.faces = s.faces(is_rim_face,:);
            
            [nrims] = find_rims(s_next);
        end
        
        if numel(nrims) > 0 % if there are nextrims
            [peaks] = find_next_rim(s, s1, Ne, nd, nrims, peaks, pki); 
        else
            if numel(next_rim) == 0
                if numel(peaks(pki).vertices) == 0
                    [~, mi] = min(nd(rim));
                    peaks(pki).apex = rim(mi);
                else
                    [~, mi] = min(nd(peaks(pki).vertices));
                    peaks(pki).apex = peaks(pki).vertices(mi);
                end
            else
                [~, mi] = min(nd(next_rim));
                peaks(pki).apex = next_rim(mi);
            end

            peaks(pki).l_center(end+1,:) = s.vertices(peaks(pki).apex,:);
            peaks(pki).height = sum(dists_points(peaks(pki).l_center));

            fprintf('peak %d %d %d %d\n',pki, parent, peaks(pki).apex, numel(peaks(pki).rims))
        end
    end
end
