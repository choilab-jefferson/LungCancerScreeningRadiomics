function log_area_param = area_distortions(mesh_orig, mesh_param)
% Normalize Mesh using ReMesh software
% Font 24

%[v_orig,f_orig] = read_obj_simple(filename_orig);
%[v_param,f_param]=read_obj_simple(filename_param);
v_orig = mesh_orig.vertices';
f_orig = mesh_orig.faces';
v_param = mesh_param.vertices';
f_param = mesh_param.faces';
area_orig = compute_per_vertex_area(v_orig,f_orig);
area_param = compute_per_vertex_area(v_param,f_param);
log_area_param = log(area_param./area_orig);
 
% area distortion histogram

%figure;
%[n,x]=hist(log_area_param, 80);
%bar(x,n./sum(n),.5,'hist');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Av = compute_per_vertex_area(V,F) 

%% Compute area around each V

% area of each face
a = V(:,F(3,:)) - V(:,F(1,:));
b = V(:,F(2,:)) - V(:,F(1,:));
ab = crossp(a',b')';
Af = sqrt(sum(ab.^2))/2;
% area of each vertex
m = size(F,2);
U = sparse( [1:m, 1:m, 1:m], [F(1,:) F(2,:) F(3,:)], [Af,Af,Af] );
Av = full(sum(U,1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = crossp(x,y)
% x and y are (m,3) dimensional
z = x;
z(:,1) = x(:,2).*y(:,3) - x(:,3).*y(:,2);
z(:,2) = x(:,3).*y(:,1) - x(:,1).*y(:,3);
z(:,3) = x(:,1).*y(:,2) - x(:,2).*y(:,1);