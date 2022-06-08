function [num_peaks,l_height,l_width,l_angle,l_nheight,l_nwidth,l_rimchange, ...
          spic_a,spic_b,spic1,spic2,spic3,spic4,hx,hn] = extract_spiculation_features(peaks, esr)
    num_peaks = numel(peaks);
    fprintf('number of peaks %d\n', num_peaks)

    if num_peaks > 0
        l_height = [peaks.height];
        l_width = [peaks.width];
        l_angle = [peaks.angle];       
        l_nheight = l_height/esr;
        l_nwidth = l_width/esr;
        
        l_height1 = [peaks.height1];
        l_spic_a = exp(-l_angle).*l_height1;
        l_spic_b = cos(l_angle).*l_height1;
        
        l_rimchange = [peaks.rimchange];
        
        spic_a = sum(l_spic_a);
        spic_b = sum(l_spic_b)/sum(l_height1);
        
        spic1 = sum([peaks.mean_ArD].*l_height)/sum(l_height);
        spic2 = sum([peaks.var_ArD].*l_height)/sum(l_height);
        spic3 = sum(l_rimchange.*l_height)/sum(l_height);
        spic4 = sum(l_rimchange(l_height>0)./l_height(l_height>0));
        %spic4 = sum([peaks.mean_ArD].*l_width)/sum(l_width);
        %spic5 = sum([peaks.var_ArD].*l_height)/sum(l_height);
        %spic6 = sum([peaks.var_ArD].*l_width)/sum(l_width);
        [hn, hx] = hist(l_nheight,0:0.1:1);
    else
        l_height = 0;
        l_width = 0;
        l_angle = 0;
        l_nheight = 0;
        l_nwidth = 0;
        
        l_rimchange = 0;
        
        spic_a = 0;
        spic_b = 0;

        spic1 = 0;
        spic2 = 0;
        spic3 = 0;
        spic4 = 0;

        hx = 0:0.1:1;
        hn = zeros(size(hx));
    end
end