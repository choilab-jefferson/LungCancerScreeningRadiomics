function [spikes] = measure_spikes(spikes,s,nd,l_ArD,attached_vertices)
    [normal_s,~] = patchnormals(s);
    num_spikes = numel(spikes);
    all_angles = compute_per_face_angles(s.vertices',s.faces');
    for pki = 1:num_spikes
        %pj = (spikes(pki).l_center(1,:)-mean(s.vertices))*(spikes(pki).l_center(2,:)-spikes(pki).l_center(1,:))';
        normal_base = mean(spikes(pki).l_center(2:end-1,:)-spikes(pki).l_center(3:end,:));
        pj = mean(normal_s([spikes(pki).apex; spikes(pki).rims{end}],:)*normal_base'> 0);
        %disp(pj)
        if pj > 0.5 % concave spike
            spikes(pki).type = 2;
        end
        
        % attached
        if numel(attached_vertices) > 0
            DD=pdist2(attached_vertices,cell2mat(spikes(pki).rims));
            spikes(pki).attachment = mean(min(DD)==0);
            if spikes(pki).attachment > 0.5
                spikes(pki).type = 3;
            end
        end
        
        if numel(spikes(pki).rims) == 1 && sum(spikes(pki).rims{1} == spikes(pki).apex) > 0 
            spikes(pki).type = 4;
            continue
        end
         
        %halfmin = nd(spikes(pki).apex)/2; % for FWHM
        halfmin = min(nd(cell2mat(spikes(pki).rims)))/2; % for FWHM
        fwhmi = 0;
        for ri = 1:numel(spikes(pki).rims)
            if mean(nd(spikes(pki).rims{ri})) < halfmin
                fwhmi = ri - 1;
                break
            end
        end
        if fwhmi > 0
            spikes(pki).fwhm_baseline = spikes(pki).rims{fwhmi};
            spikes(pki).fwhm_baseline_curve = s.vertices(spikes(pki).fwhm_baseline([1:end,1]),:);
            
            spikes(pki).fwhm_height = sum(dists_points(spikes(pki).l_center(fwhmi:end,:)));
        else
            spikes(pki).fwhm_baseline = spikes(pki).baseline;
            spikes(pki).fwhm_baseline_curve = spikes(pki).baseline_curve;
            
            spikes(pki).fwhm_height = spikes(pki).height/2;
        end
        
        
        %% height
        spikes(pki).height = spikes(pki).fwhm_height * 2;
        spikes(pki).height1 = sqrt(sum((spikes(pki).l_center(end,:)-spikes(pki).l_center(1,:)).^2,2));
        
        %% angle ditribution of the faces 
        sel_face = sum(ismember(s.faces,spikes(pki).vertices),2)>0;
        angles = min(all_angles(:,sel_face));
        spikes(pki).angles = angles(:)';
        if size(spikes(pki).fwhm_baseline) < 3
            spikes(pki).type = 4;
            spikes(pki).width = 0;
            spikes(pki).angle = 0;
            continue; 
        end
        
        %% width
        w = max(pdist(s.vertices(spikes(pki).fwhm_baseline,:)));
        spikes(pki).width=w;
        
        %% angle of the spike
        nb = numel(spikes(pki).fwhm_baseline);
%        st.vertices = [s.vertices(spikes(pki).fwhm_baseline,:);mean(s.vertices(spikes(pki).fwhm_baseline,:),1);s.vertices(spikes(pki).apex,:)];
%        st.faces = [(1:nb)',repmat(nb+1,nb,1),repmat(nb+2,nb,1)];
%        fanglet=zeros(nb,3);
%         %% largest angle calculation
%         for fi = 1:nb
%             fanglet(fi,:) = triangle_angles(st.vertices(st.faces(fi,:),:));
%         end
%         fanglet(isnan(fanglet)) = [];
%         try
%             %spikes(pki).angle=mean(fanglet(:,3)*2);
%             spikes(pki).angle=max(fanglet(:,3))+min(fanglet(:,3));
%         catch
%             spikes(pki).angle=0;
%         end
%         
        %% solid angle calculation
        st.vertices = [s.vertices(spikes(pki).fwhm_baseline,:);mean(s.vertices(spikes(pki).fwhm_baseline,:),1);s.vertices(spikes(pki).apex,:)];
        st.faces = [(1:nb)',([2:nb,1])',repmat(nb+1,nb,1)];
        farealet=zeros(nb,1);
        for fi = 1:nb
            farealet(fi) = triangle_area(st.vertices(st.faces(fi,:),:));
        end
        farealet(isnan(farealet)) = [];
        try
            spikes(pki).angle=sum(farealet)/(spikes(pki).height)^2;
        catch
            spikes(pki).angle=0;
        end
        
        if spikes(pki).type == 0 && spikes(pki).angle > 0.65 % lobulaiton
            spikes(pki).type = 1;
        end
    end
    
    
    %% statistics of area distortion
    for pki = 1:num_spikes
        V = unique([cell2mat(spikes(pki).rims);spikes(pki).vertices;spikes(pki).apex]);
        spikes(pki).mean_ArD = mean(l_ArD(V));
        spikes(pki).var_ArD = var(l_ArD(V));
        spikes(pki).rimchange = sum(dists_points(s.vertices(spikes(pki).baseline([1:end,1]),:)))-sum(dists_points(s.vertices(spikes(pki).rims{end}([1:end,1]),:)));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = compute_per_face_angles(X,F)
    m = size(F,2);
    A = zeros(3,m);

    for i=1:3
       i2 = mod(i  ,3)+1;
       i3 = mod(i+1,3)+1;
       u = X(:,F(i2,:)) - X(:,F(i,:));
       v = X(:,F(i3,:)) - X(:,F(i,:));
       % normalize the vectors   
       u = u ./ repmat( sqrt(sum(u.^2,1)), [3 1] );
       v = v ./ repmat( sqrt(sum(v.^2,1)), [3 1] );
       % compute angles
       A(i,:) = acos( sum( u.*v, 1 ) );
    %    alpha = 1 ./ tan( acos(sum(u.*v,1)) );
    %    alpha = max(alpha, 1e-2); % avoid degeneracy
    end
end
