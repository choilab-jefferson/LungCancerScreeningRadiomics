function [n x x2] = mean_proj(input,input2,weight,range,range2)
%     additional = sum(input > range(end)) > 0;
%     sumw = sum(weight);
%     
%     tick = (range(end) - range(1))/(size(range,2)-1);
%     
%     for i = 1:size(range,2) + additional
%         if i == 1
%             minr = -Inf;
%             maxr = range(i);
%             x(i) = range(i) - tick/2;
%         elseif i == size(range,2) + 1
%             minr = range(i-1);
%             maxr = Inf;
%             x(i) = range(i-1) + tick/2;
%         else
%             minr = range(i-1);
%             maxr = range(i);
%             x(i) = range(i-1) + tick/2;
%         end
%         
%         n(i) = sum((input > minr & input <= maxr) .* weight);
%     end

    edges = range(:)';
    binwidth = diff(edges);
    % Set edges for call to bar below.
    edges = [-Inf, edges(1:end-1)+binwidth/2, Inf]; 
    
    edges2 = range2(:)';
    binwidth2 = diff(edges2);
    % Set edges for call to bar below.
    edges2 = [-Inf, edges2(1:end-1)+binwidth2/2, Inf]; 
    
    n = zeros(size(edges,2)-1, size(edges2,2)-1);
    for i = 1:size(edges,2) - 1
        minr = edges(i);
        maxr = edges(i+1);
        x(i) = range(i);
        for j = 1:size(edges2,2) - 1
            minr2 = edges2(j);
            maxr2 = edges2(j+1);
            x2(j) = range2(j);
            
            t = weight((input > minr & input <= maxr)&(input2 > minr2 & input2 <= maxr2));
            if sum(abs(t(:))) > 0
                n(i,j) = mean(t(abs(t)>0));
            end
        end
    end
end
