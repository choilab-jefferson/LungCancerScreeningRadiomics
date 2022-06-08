function [s1,colors,normals] = sphereical_mapping(s)
    n = size(s.vertices,1);
	m = size(s.faces,1);
    ce = mean(s.vertices);
    v = s.vertices - repmat(ce,n,1);
    magv = sqrt(sum(v.^2,2));
    mdist = mean(magv);
    nv = v./repmat(magv,1,3);

    %% Compute the weights. The weights should be positive for the method to work.
    W = sparse(n,n);
    for i=1:3
        i1 = mod(i-1,3)+1;
        i2 = mod(i  ,3)+1;
        i3 = mod(i+1,3)+1;
        pp = s.vertices(s.faces(:,i2),:) - s.vertices(s.faces(:,i1),:);
        qq = s.vertices(s.faces(:,i3),:) - s.vertices(s.faces(:,i1),:);
        % normalize the vectors
        pp = pp ./ repmat( sqrt(sum(pp.^2,2)), [1 3] );
        qq = qq ./ repmat( sqrt(sum(qq.^2,2)), [1 3] );
        % compute angles
        ang = acos(sum(pp.*qq,2));
        a = max(1 ./ tan(ang),1e-1); % this is *very* important
        W = W + sparse(s.faces(:,i2),s.faces(:,i3), a, n, n );
        W = W + sparse(s.faces(:,i3),s.faces(:,i2), a, n, n );
    end
    %W(abs(W)>0) = 1;

    %% Compute the normalized weight matrix tW such that its rows sums to 1.
    d = full( sum(W,1) );
    D = spdiags(d(:), 0, n,n);
    iD = spdiags(d(:).^(-1), 0, n,n);
    tW = iD * W;


    %% Compute an initial mapping on the sphere. This simply a radial projection.
    s1 = s;
    s1.vertices = nv;


    %% Check which faces have the correct orientation.
    % normal to faces
    [~,normalf] = patchnormals(s1);
    % center of faces
    C = squeeze(mean(reshape(s1.vertices(s1.faces,:),[m 3 3]), 2));
    % inner product
    I = sum(C.*normalf,2);
    
    if sum(I(:)<0)/m > 0.5
        s.faces = s.faces(:,[2 1 3]);
        s1.faces = s.faces;
        normalf = -normalf;
        I = sum(C.*normalf,2);
    end

    disp(['Ratio of inverted triangles:' num2str(sum(I(:)<0)/m*100, 3) '%']);

    %% Perform smoothing and projection.
    niter = 1000;
    ndisp = round([0 0.1 0.3 1]*niter); ndisp = max(ndisp,1);
    k = 1;

    Edir = zeros(niter,1);
    for i=1:niter
        % smooth 
        s1.vertices = tW*s1.vertices;
        % project
        s1.vertices = s1.vertices ./ repmat( sqrt(sum(s1.vertices.^2,2)), [1 3] );
        
        
        %Dirichlet energy
        E = 0;
        for j=1:3
            j1 = mod(j,3)+1;
            % directed edge
            u = s1.vertices(s1.faces(:,j),:) - s1.vertices(s1.faces(:,j1),:);
            % norm squared
            u = sum(u.^2,2);
            % weights between the vertices
            w = W(s1.faces(:,j) + (s1.faces(:,j1)-1)*n);
            E = sum( w.*u );
        end
        Edir(i) = E;
        
         % compute inverted triangles
        [normals,normalf] = patchnormals(s1);
        C = squeeze(mean(reshape(s1.vertices(s1.faces,:),[m 3 3]), 2));
        I = sum(C.*normalf,2);
        ninvert = sum(I<0);
        
        if ninvert==0 || i > 1000 && Edir(i-1)-Edir(i) < 1e-16
            %disp([i, Edir(i),  Edir(i-1)])
            break
        end
    end


    disp(['Ratio of inverted triangles:' num2str(ninvert/m*100, 3) '%: err:' num2str(Edir(1)) '->' num2str(Edir(end))]);
    
    s1.vertices = s1.vertices*mdist + repmat(ce,n,1);
    ce1 = mean(s1.vertices);
    v1 = s1.vertices-repmat(ce1,n,1);
    magv1 = sqrt(sum(v1.^2,2));
    nv1 = v1./repmat(magv1,1,3);
    
   
    colors = (1-normals)/2; % recaculate vetex color
    
