function path_ij = find_shortestpath(s, i, j)
    global geodesic_library;
    geodesic_library = 'geodesic_release';
    mesh = geodesic_new_mesh(s.vertices,s.faces);         %initilize new mesh
    algorithm = geodesic_new_algorithm(mesh, 'dijkstra');      %initialize new geodesic algorithm
    source_points = {geodesic_create_surface_point('vertex',i,s.vertices(i,:))};
    geodesic_propagate(algorithm, source_points);   %propagation stage of the algorithm (the most time-consuming)
    
    destination = geodesic_create_surface_point('vertex',j,s.vertices(j,:));
    path = geodesic_trace_back(algorithm, destination);     %find a shortest path from source to destination
    
    path_ij = zeros(length(path),1);
    for idx=1:length(path)
        path_ij(idx) = path{idx}.id;
    end;
    path_ij = flip(path_ij);
end