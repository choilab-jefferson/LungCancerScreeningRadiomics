function [spikes] = detect_spikes(s, s1, Ne, nd, spikes)
%     patch(s, 'FaceVertexCData',nd,'FaceAlpha',0.4,'FaceColor','interp', 'EdgeColor', 'none');
%     colormap('jet'), caxis([-1 1])
%     view(-37.5,30);
%     axis equal fill
%     hold on

    global geodesic_library;                
    geodesic_library = 'geodesic_release';      %"release" is faster and "debug" does additional checks
    
    for i = 1:numel(spikes)
        spikes = find_next_rim(s, s1, Ne, nd, {spikes(i).baseline}, spikes, i);
    end
    idx = [spikes.apex]>0;
    spikes = spikes(idx);
    
    geodesic_delete;
end

function [spikes] = find_next_rim(s, s1, Ne, nd, rims, spikes, parent)
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
        
        pki = numel(spikes)+1;
        spikes(pki) = spikes(parent);
        spikes(pki).parent = parent;
        
        % when there are multiple rims and its perimeter is less than 60% of entire perimeter
        % replace baseline
        if 1 && nr > 1 && p_rims(rid)/sum(p_rims) < 0.6 && dists_points([spikes(pki).center; center]) > 2*spikes(parent).height/numel(spikes(parent).rims)
            spikes(pki).baseline = spikes(parent).rims{end};
            spikes(pki).spike_vertices = spikes(parent).vertices;
            spikes(pki).baseline_curve = s.vertices(rim([1:end,1]),:);
            spikes(pki).rims = {rim};
            spikes(pki).l_center = center;
            spikes(pki).l_normal = normal;
            fprintf('baseline %d %d %d %0.2f\n',pki, parent, nr, p_rims(rid))
        else % add a rim into the spike
            spikes(pki).rims = [spikes(pki).rims;rim];
            spikes(pki).l_center = [spikes(pki).l_center;center];
            spikes(pki).l_normal = [spikes(pki).l_normal;normal];
        end

        % update spike measures
        spikes(pki).height = sum(dists_points(spikes(pki).l_center));
        spikes(pki).center = center;
        spikes(pki).normal = normal;
                
        
        %% find inside vertices
        v_rim1 = s1.vertices(rim,:);
        v_spike1 = s1.vertices(spikes(pki).vertices,:);
        [~,normal1] = calcaulte_center_normal(v_rim1);
        xy_rim = RodriguesRotation(v_rim1,normal1,[0, 0, 1]);
        xy_spike = RodriguesRotation(v_spike1,normal1,[0, 0, 1]);
        dist2rim = min(pdist2(v_rim1,v_spike1));
        in_vertices = inpolygon(xy_spike(:,1),xy_spike(:,2),xy_rim(:,1),xy_rim(:,2));
        in_vertices1 = (dist2rim>0&dist2rim<max(pdist(v_rim1)))'; % not in the rim
                
        spikes(pki).vertices = spikes(pki).vertices(in_vertices&in_vertices1); % select inside vertices

        medndr = median(nd(rim));
        mndr = mean(nd(rim));
        %next_rim = spikes(pki).vertices(nd(spikes(pki).vertices)<min(medndr,mndr));
        next_rim = spikes(pki).vertices(nd(spikes(pki).vertices)<min(nd(rim)));

%         v_inside = s.vertices(spikes(pki).vertices,:);
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
            [spikes] = find_next_rim(s, s1, Ne, nd, nrims, spikes, pki); 
        else
            if numel(next_rim) == 0
                if numel(spikes(pki).vertices) == 0
                    [~, mi] = min(nd(rim));
                    spikes(pki).apex = rim(mi);
                else
                    [~, mi] = min(nd(spikes(pki).vertices));
                    spikes(pki).apex = spikes(pki).vertices(mi);
                end
            else
                [~, mi] = min(nd(next_rim));
                spikes(pki).apex = next_rim(mi);
            end

            spikes(pki).l_center(end+1,:) = s.vertices(spikes(pki).apex,:);
            spikes(pki).height = sum(dists_points(spikes(pki).l_center));

            fprintf('spike %d %d %d %d\n',pki, parent, spikes(pki).apex, numel(spikes(pki).rims))
        end
    end
end
