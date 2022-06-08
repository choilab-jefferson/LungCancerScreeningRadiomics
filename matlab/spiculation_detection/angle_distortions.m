function log_angle_param = angle_distortions(mesh_orig, mesh_param)
% Normalize Mesh using ReMesh software
% Font 24

%[v_orig,f_orig] = read_obj_simple(filename_orig);
%[v_param,f_param]=read_obj_simple(filename_param);
v_orig = mesh_orig.vertices';
f_orig = mesh_orig.faces';
v_param = mesh_param.vertices';
f_param = mesh_param.faces';
angle_orig = compute_per_face_min_angle(v_orig,f_orig);
angle_param = compute_per_face_min_angle(v_param,f_param);
log_angle_param = log(angle_param./angle_orig);
 
% % area distortion histogram
% 
% figure;
% [n,x]=hist(log_angle_param, 80);
% bar(x,n./sum(n),.5,'hist');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = compute_per_face_min_angle(X,F)

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
A = min(A);
