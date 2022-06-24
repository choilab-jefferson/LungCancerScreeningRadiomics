function [spikes] = detect_baselines(s, nd)
    % initialize spikes
    spikes = struct('baseline',{},'vertices',{},'rims',{},'l_center',{},'l_normal',{},'apex',{},'height',{},'width',{},'angle',{});
    
    s_base = s;
    n_negative = sum(nd(s.faces) <= 0,2); % count negative area distortion vertices for each face
    s_base.faces = s.faces(n_negative>0,:); % if any can be a part of spike
    baseline_f = find(n_negative>0);
   
    %% labeling vertices of baseline faces
    baseline_v = unique(s_base.faces(:));
    Ne = vertex_neighbours(s_base);
    n = size(baseline_v,1);
    label = zeros(n,1);
    ci = 1;
    rims = {};
    for bi = 1:n
        if label(bi)~=0
            continue
        else
            [label] = labeling(baseline_v,Ne,ci,bi,label);
            % idx = baseline_v(label==ci);
            % plot3(s_base.vertices(idx,1),s_base.vertices(idx,2),s_base.vertices(idx,3),'*'), hold on
            ci = ci + 1;
        end
    end
    
    %% convert vertax label to face labels
    m = size(s_base.faces,1);
    flabel = zeros(m,1);
    for fi = 1:m
        flabel(fi) = max(label((baseline_v == s_base.faces(fi,1))+(baseline_v == s_base.faces(fi,2))+(baseline_v == s_base.faces(fi,3))>0));
    end
    
    %% find baselines
    pki = 1;
    for ci = 1:max(label)
        s_base_ci = s_base;
        s_base_ci.faces = s_base.faces(flabel == ci,:);
        [rims] = find_rims(s_base_ci);
        nr = numel(rims);
        for ri = 1:nr
            rim = rims{ri};
            spikes(pki).type = 0;
            spikes(pki).baseline = rim;
            spikes(pki).baseline_curve = s.vertices(rim([1:end,1]),:);            
            spikes(pki).apex = 0;
            spikes(pki).rims = {rim};
            spikes(pki).vertices = baseline_v(label == ci);
            spikes(pki).spike_vertices = spikes(pki).vertices;
            spikes(pki).faces = baseline_f(flabel == ci);
            [center,normal] = calcaulte_center_normal(s.vertices(rim,:));
            spikes(pki).l_center = [spikes(pki).l_center;center];
            spikes(pki).l_normal = [spikes(pki).l_normal;normal];
            spikes(pki).center = center;
            spikes(pki).normal = normal;
            spikes(pki).angle = 0;
            spikes(pki).height = 0;
            spikes(pki).height1 = 0;
            spikes(pki).width = 0;
            
            pki=pki+1;
            % plot3(s.vertices(rim([1:end,1]),1),s.vertices(rim([1:end,1]),2),s.vertices(rim([1:end,1]),3),'-'), hold on
        end
    end
end

function [label] = labeling(V,Ne,c,bi,label)
    label(bi) = c;
    Neb = Ne{V(bi)};
    while(numel(Neb)>0)
        ne = Neb(1);
        for nbi = find(ne == V)
           if label(nbi)==0
                % disp([bi, nbi, c])
                label(nbi) = c;
                Nenb = Ne{V(nbi)};
                Neb(end+1:end+numel(Nenb)) = Nenb;
            end
        end
        Neb(1) = []; % delete
    end
end


