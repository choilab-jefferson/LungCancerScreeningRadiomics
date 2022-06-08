function [peaks] = measure_peaks(peaks,s,nd,l_ArD,attached_vertices)
    [normal_s,~] = patchnormals(s);
    num_peaks = numel(peaks);
    all_angles = compute_per_face_angles(s.vertices',s.faces');
    for pki = 1:num_peaks
        %pj = (peaks(pki).l_center(1,:)-mean(s.vertices))*(peaks(pki).l_center(2,:)-peaks(pki).l_center(1,:))';
        normal_base = mean(peaks(pki).l_center(2:end-1,:)-peaks(pki).l_center(3:end,:));
        pj = mean(normal_s([peaks(pki).apex; peaks(pki).rims{end}],:)*normal_base'> 0);
        %disp(pj)
        if pj > 0.5 % concave peak
            peaks(pki).type = 2;
        end
        
        % attached
        if numel(attached_vertices) > 0
            DD=pdist2(attached_vertices,cell2mat(peaks(pki).rims));
            peaks(pki).attachment = mean(min(DD)==0);
            if peaks(pki).attachment > 0.5
                peaks(pki).type = 3;
            end
        end
        
        if numel(peaks(pki).rims) == 1 && sum(peaks(pki).rims{1} == peaks(pki).apex) > 0 
            peaks(pki).type = 4;
            continue
        end
         
        %halfmin = nd(peaks(pki).apex)/2; % for FWHM
        halfmin = min(nd(cell2mat(peaks(pki).rims)))/2; % for FWHM
        fwhmi = 0;
        for ri = 1:numel(peaks(pki).rims)
            if mean(nd(peaks(pki).rims{ri})) < halfmin
                fwhmi = ri - 1;
                break
            end
        end
        if fwhmi > 0
            peaks(pki).fwhm_baseline = peaks(pki).rims{fwhmi};
            peaks(pki).fwhm_baseline_curve = s.vertices(peaks(pki).fwhm_baseline([1:end,1]),:);
            
            peaks(pki).fwhm_height = sum(dists_points(peaks(pki).l_center(fwhmi:end,:)));
        else
            peaks(pki).fwhm_baseline = peaks(pki).baseline;
            peaks(pki).fwhm_baseline_curve = peaks(pki).baseline_curve;
            
            peaks(pki).fwhm_height = peaks(pki).height/2;
        end
        
        
        %% height
        peaks(pki).height = peaks(pki).fwhm_height * 2;
        peaks(pki).height1 = sqrt(sum((peaks(pki).l_center(end,:)-peaks(pki).l_center(1,:)).^2,2));
        
        %% angle ditribution of the faces 
        sel_face = sum(ismember(s.faces,peaks(pki).vertices),2)>0;
        angles = min(all_angles(:,sel_face));
        peaks(pki).angles = angles(:)';
        if size(peaks(pki).fwhm_baseline) < 3
            peaks(pki).type = 4;
            peaks(pki).width = 0;
            peaks(pki).angle = 0;
            continue; 
        end
        
        %% width
        w = max(pdist(s.vertices(peaks(pki).fwhm_baseline,:)));
        peaks(pki).width=w;
        
        %% angle of the peak
        nb = numel(peaks(pki).fwhm_baseline);
%        st.vertices = [s.vertices(peaks(pki).fwhm_baseline,:);mean(s.vertices(peaks(pki).fwhm_baseline,:),1);s.vertices(peaks(pki).apex,:)];
%        st.faces = [(1:nb)',repmat(nb+1,nb,1),repmat(nb+2,nb,1)];
%        fanglet=zeros(nb,3);
%         %% largest angle calculation
%         for fi = 1:nb
%             fanglet(fi,:) = triangle_angles(st.vertices(st.faces(fi,:),:));
%         end
%         fanglet(isnan(fanglet)) = [];
%         try
%             %peaks(pki).angle=mean(fanglet(:,3)*2);
%             peaks(pki).angle=max(fanglet(:,3))+min(fanglet(:,3));
%         catch
%             peaks(pki).angle=0;
%         end
%         
        %% solid angle calculation
        st.vertices = [s.vertices(peaks(pki).fwhm_baseline,:);mean(s.vertices(peaks(pki).fwhm_baseline,:),1);s.vertices(peaks(pki).apex,:)];
        st.faces = [(1:nb)',([2:nb,1])',repmat(nb+1,nb,1)];
        farealet=zeros(nb,1);
        for fi = 1:nb
            farealet(fi) = triangle_area(st.vertices(st.faces(fi,:),:));
        end
        farealet(isnan(farealet)) = [];
        try
            peaks(pki).angle=sum(farealet)/(peaks(pki).height)^2;
        catch
            peaks(pki).angle=0;
        end
        
        if peaks(pki).type == 0 && peaks(pki).angle > 0.65 % lobulaiton
            peaks(pki).type = 1;
        end
    end
    
    
    %% statistics of area distortion
    for pki = 1:num_peaks
        V = unique([cell2mat(peaks(pki).rims);peaks(pki).vertices;peaks(pki).apex]);
        peaks(pki).mean_ArD = mean(l_ArD(V));
        peaks(pki).var_ArD = var(l_ArD(V));
        peaks(pki).rimchange = sum(dists_points(s.vertices(peaks(pki).baseline([1:end,1]),:)))-sum(dists_points(s.vertices(peaks(pki).rims{end}([1:end,1]),:)));
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
