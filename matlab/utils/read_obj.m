function [vertex,face] = read_obj(input_file)

%   [vertex,face,color,normal] = read_obj(filename);
%
%   faces    : list of facesangle elements
%   vertex  : node vertexinatates
%   normal : normal vector list

fid = fopen(input_file);
if fid<0
    error(['Cannot open ' input_file '.']);
end

frewind(fid);
normal = [];
vertex = [];
face = [];
while 1
    s = fgetl(fid);
    if ~ischar(s), 
        break;
    end
    if ~isempty(s) && strcmp(s(1), 'f')
        % face
        s = strrep(s,'//','/');
        face(:,end+1) = sscanf(s(3:end), '%d/%d %d/%d %d/%d');
    end
    if ~isempty(s) && strcmp(s(1), 'v')
        % vertex
        if ~isempty(s) && strcmp(s(2), 'n')
            normal(:,end+1) = sscanf(s(3:end), '%f %f %f');
        elseif ~isempty(s) && strcmp(s(2), 't')
            normal(:,end+1) = sscanf(s(3:end), '%f %f');
        else
            vertex(:,end+1) = sscanf(s(3:end), '%f %f %f');
        end
    end
end
fclose(fid);

face = [face(1,:); face(3,:); face(5,:)]';
vertex = vertex(1:3,:)';