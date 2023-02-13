function spherical_param(sph_map_filename, obj_filename)
    s = struct();
    s1 = s;
    [s.vertices,s.faces] = readOBJ(obj_filename);
    [~,name,~] = fileparts(obj_filename);
    tempheader = tempname;
    mkdir(tempheader);
    temp_ply = [tempheader '/' name '.ply'];
    disp(temp_ply)
    out_header = [tempheader '/out'];
    writePLY(temp_ply, s.vertices, s.faces);
    disp('Spherical parameterization...')
    if isfile('/usr/bin/docker')
        ConformalizedMCF = ['docker run --user 1007 -v ' tempheader ':' tempheader ' wookjinchoi/conformalized_mcf:latest ConformalizedMCF'];
    else 
        ConformalizedMCF = ['ConformalizedMCF'];
    end
    system([ConformalizedMCF ' --in ' temp_ply ' --outHeader ' out_header ' --steps 3000 --threads 1']);
    [s1.vertices,s1.faces] = readPLY([out_header '.3000.ply']);
    writeOBJ(sph_map_filename, s1.vertices, s1.faces);
    delete([tempheader '/*']);
    [status, message, messageid] = rmdir(tempheader);
    if status == 0
        disp(message)
        disp(messageid)
    end
end
