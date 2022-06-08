function [center,normal] = calcaulte_center_normal(v_baseline)
    center = mean(v_baseline,1);
    nb = size(v_baseline,1);
    v_baseline_to_origin = v_baseline - ones(nb, 1)*center;
    [~, ~, V]= svd(v_baseline_to_origin);
    normal = V(:, 3)';
end