function [vertex,face,normal] = read_obj(filename)

%   [vertex,face,normal,texture] = read_obj(filename);
%
%   faces    : list of facesangle elements
%   vertex  : node vertexinatates
%   normal : normal vector list
%
%   Copyright (c) 2008 Gabriel Peyre

fid = fopen(filename);
if fid<0
    error(['Cannot open ' filename '.']);
end

frewind(fid);
a = fscanf(fid,'%c',1);
if strcmp(a, 'P')
    % This is the montreal neurological institute (MNI) specific ASCII facesangular mesh data structure.
    % For FreeSurfer software, a slightly different data input coding is
    % needed. It will be provided upon request.
    fscanf(fid,'%f',5);
    n_points=fscanf(fid,'%i',1);
    vertex=fscanf(fid,'%f',[3,n_points]);
    normal=fscanf(fid,'%f',[3,n_points]);
    n_faces=fscanf(fid,'%i',1);
    fscanf(fid,'%i',5+n_faces);
    face=fscanf(fid,'%i',[3,n_faces])'+1;
    fclose(fid);
    return;
end

frewind(fid);
vertex = [];
face = [];
normal = [];
texture = [];
while 1
    s = fgetl(fid);
    if ~ischar(s), 
        break;
    end
    if ~isempty(s) && strcmp(s(1), 'f')
        % face
        face(:,end+1) = sscanf(s(3:end), '%d//%d %d//%d %d//%d');
    end
    if ~isempty(s) && strcmp(s(1), 'v')
        % vertex
        if ~isempty(s) && strcmp(s(2), 'n')
            normal(:,end+1) = sscanf(s(3:end), '%f %f %f');
        else
            vertex(:,end+1) = sscanf(s(3:end), '%f %f %f');
        end
        %vertex(:,end+1) = sscanf(s(3:end), '%f %f %f %f %f %f');
        %vertex(:,end+1) = sscanf(s(3:end), '%f %f %f');
    end
end
face = [face(1,:); face(3,:); face(5,:)];
fclose(fid);