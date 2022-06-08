function Mapnew=makenewscalarmap(Map,HT_index,HT_values)
% Combine the edge middle points and old vertex points to faces.
% (4 Split method)
% Make output V
Mapout=zeros(size(HT_index,1)*4,1);
Mapout(1:size(HT_index,1))=Map;

max_Vindex = 0;
for i=1:length(HT_index)
    HT_indexi = HT_index{i};
    for j=1:length(HT_indexi)
        P0 = Map(i);
        P1 = Map(HT_indexi(j));
        Pnew = (P0+P1)/2;
        Vindex = HT_values{i}(j);
        %if Mapout(Vindex) == 0
        Mapout(Vindex,:)=Pnew;
        %else
            %disp([Mapout(Vindex,:), Pnew])
        %end
        if Vindex > max_Vindex
            max_Vindex = Vindex;
        end
    end
end

Mapnew=Mapout(1:max_Vindex);
