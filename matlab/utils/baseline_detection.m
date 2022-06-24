function [spikes,t]=baseline_detection(spikes, s, nd, Ne)
    n_spikes = numel(spikes);
    t = zeros(size(nd));
    nd1 = nd;
    for vi = 1:numel(nd)
        nd1(vi) = mean(nd([Ne{vi};vi]));
    end
    for pki = 1:n_spikes % check negibourhoond
        count = 1;
        t1 = zeros(size(nd));
        ppi = spikes(pki).apex;
        %if t(ppi)>0, continue; end
        t1(ppi) = 1;
        ne = cell2mat(Ne(ppi));
        while numel(ne)>0
            count = count + 1;
            t1(ne)=1;
            bt1 = t1==3;

            nne = unique(cell2mat(Ne(ne)));
            if sum(bt1)/numel(nne)>2
                break
            end
            nne = nne(t1(nne)==0);

            for nei = ne'
                nnei = cell2mat(Ne(nei));
                nnei = nnei(t1(nnei)==0);
                pv = nd1(nei);
                nv = nd1(nnei);
                bt = nv > 0;
                %if pv>nd1(ppi)+0.2 && mean(pv-nv>0) > 0.5
                %    t1(nnei) = 2;
                %end
                if sum(bt)
                    t1(nei) = 2;
                    t1(nnei) = 3;
                    %ndisp([pki, count, pv, nv', sum(bt)])
                end
            end

            ne = nne(t1(nne)==0);
        end

        %spikes(pki).baseline = [s.vertices(t1==3,:);s.vertices(t1==2,:)];
        spikes(pki).baseline = find(t1==2);
        spikes(pki).vertices = find(t1>0&t1<3)';
        %spike_v{pki} = nd(t1>0);

        t(ppi) = 4;
        t=bitor(t,t1);
    end