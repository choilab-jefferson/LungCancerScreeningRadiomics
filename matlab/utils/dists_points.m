function dists = dists_points(points)
    dists = sqrt(sum((points(1:end-1,:)-points(2:end,:)).^2,2));
end