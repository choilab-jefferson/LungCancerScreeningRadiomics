function [fig] = plot_contours_and_attachment1(c_img, seg_img, cip_img, gc_img, att_img)
    fig = figure('InvertHardcopy','off','Color',[1 1 1],'visible','off','PaperUnits', 'inches','PaperPosition', [0 0 4 4]);
    I = gray2rgb(uint8(double(c_img).*(~att_img)),gray)+gray2rgb(uint8(double(c_img)/4*3.*(att_img==255)),gray);
    I =  I + gray2rgb(att_img,[0,0,0;0,1,0])/4*1;
    imshow(I); hold on; 
    s_b = bwboundaries(seg_img); c_b = bwboundaries(cip_img); g_b = bwboundaries(gc_img);
    for k = 1:length(g_b)
        p = plot(g_b{k}(:,2),g_b{k}(:,1),'b','LineWidth',1.5);
        p.Color(4) = 0.3;
    end
    for k = 1:length(c_b)
        p = plot(c_b{k}(:,2),c_b{k}(:,1),'r','LineWidth',1.5);
        p.Color(4) = 0.3;
    end
    for k = 1:length(s_b)
        p = plot(s_b{k}(:,2),s_b{k}(:,1),'y','LineWidth',3);
        p.Color(4) = 0.7;
    end
    set(gca,'Position',[0 0 1 1],'units','normalized');
end