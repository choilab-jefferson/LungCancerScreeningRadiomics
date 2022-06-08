function [rims] = find_rims(s)
    F = s.faces;
    V = s.vertices;
    % Find all edges in mesh, note internal edges are repeated
    E = sort([F(:,1) F(:,2); F(:,2) F(:,3); F(:,3) F(:,1)]')';
    % determine uniqueness of edges
    [u,m,n] = unique(E,'rows');
    % determine counts for each unique edge
    counts = accumarray(n(:), 1);
    % extract edges that only occurred once
    O = u(counts==1,:);
    
    n_edges = size(O,1);
    label = zeros(n_edges,1);
    ci = 1;
    rims = {};
    while(n_edges > sum(label>0))
        rim = [];
        [label,rim] = labeling(O,ci,label,rim);
        if numel(rim) < 4
            continue
        end
        if 0
            rims = [rims;rim];
            ci = ci + 1;
        else % rim cut
            nrims = cut_rims(rim, s, rim);
            rims = [rims;nrims];
            ci = ci + numel(nrims);
        end
    end
end

function [label, rim] = labeling(O,c,label,rim)
    u = O(~label,:);
    n = u(1,:);
    while numel(n) > 0
        e = n(1,:);
        rim = [rim; e(2)];
        label(O(:,1)==e(1)&O(:,2)==e(2)|O(:,2)==e(1)&O(:,1)==e(2)) = c;
        
        u = O(~label,:);
        n = u(u(:,1) == e(2),:);
        if numel(n) == 0
            n2 = u(u(:,2) == e(2),:);
            n = n2(:,[2 1]);
        end
    end
end

function [D1,D2,D3] = rim_distances(rim, s)
    D1 = pdist2(s.vertices(rim,:),s.vertices(rim,:)); % minimum distance of vertices pairs
    n = size(rim,1);
    D2 = nan(n);
    D3 = nan(n);
    for i = 1:n
        for j = i+2:(n-floor(2/i))
            c_cwh = rim(i:j);
            c_ccwh = rim([j:end,1:i]);
            
            D2(i,j) = sum(sqrt(sum((s.vertices(c_cwh(2:end),:)-s.vertices(c_cwh(1:end-1),:)).^2,2))); % clock wise arc length
            D3(i,j) = sum(sqrt(sum((s.vertices(c_ccwh(2:end),:)-s.vertices(c_ccwh(1:end-1),:)).^2,2))); % counter clock wise arc length
        end
    end
    D1(isnan(D2))=nan;
end

function rims = cut_rims(rim, s, rim1)
    rims = {};
    [D1,D2,D3] = rim_distances(rim1, s);
    D4 = D1./min(D2,D3); % shortest distance / minor arc length
    D5 = D1./max(D2,D3); % shortest distance / major arc length
    [mn,k] = min(D4(:)); % min distance pair of minor arc ratio
    io = mod(k,numel(rim1)); % start index
    jo = floor(k/numel(rim1))+1; % end index
    mx = D5(io,jo); % distance of major arc length with min distance minor arc
    rim_rh1 = rim1([1:io,jo:end]); % right arc
    rim_lh1 = rim1(io:jo); % left arc
    i = find(rim == rim1(io),1,'first');
    j = find(rim == rim1(jo),1,'last');
    
    if mn < 0.3 && mx < 0.1 && numel(rim_rh1)>3 && numel(rim_lh1)>3
        try
            path_ij = find_shortestpath(s,rim(i),rim(j));
            rim_rh = [rim(1:i-1);path_ij;rim(j+1:end)];
            rim_lh = [rim(i+1:j-1);flip(path_ij)];
            if numel(rim_rh1)>3
                rims_rh = cut_rims(rim_rh,s,rim_rh1);
                rims = [rims;rims_rh];
            end
            if numel(rim_lh1)>3
                rims_lh = cut_rims(rim_lh,s,rim_lh1);
                rims = [rims;rims_lh];
            end
            fprintf('rims %d %d %d\n',numel(rims),numel(rim_rh1),numel(rim_lh1))
        catch
            rims = {rim};
        end
    else
        rims = {rim};
    end
end