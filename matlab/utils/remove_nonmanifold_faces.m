function [SV,SF] = remove_nonmanifold_faces(V,F)        
    SV = V;
    nme = nonmanifold_edges(F);
    remove = logical(size(F,1));
    for i = 1:size(F,1)
        remove(i) =  sum(sum(F(i,:) == unique(nme(:)))) > 0;
    end
    SF = F(~remove,:);
    [SV,SF] = remove_degenerate_faces(SV,SF);
    HF=fill_holes(SV,SF);
    SF = [HF;SF];
end