%% main function
clc;
clear;
close all

set_environment; % call environment setting

%% spiculation toolbox
addpath(genpath('utils'))
addpath(genpath('spiculation_detection'))

%% set values
experiment_set = 'nodule-lidc';

%% directory paths
experiment_path = [data_path '/' experiment_set];
output_experiment_path=['../' output_path '/' experiment_set];
if ~isdir(output_experiment_path); mkdir(output_experiment_path); end
if ~isdir([output_experiment_path '/spiculation']); mkdir([output_experiment_path '/spiculation']); end
if ~isdir([output_experiment_path '/parameters']); mkdir([output_experiment_path '/parameters']); end

obj_path = [experiment_path '/objs'];
sph_map_path = [experiment_path '/spherical_obj'];
dir_obj = dir([obj_path '/*.obj']);
pid_list = {dir_obj.name};

features = table();
all_spikes = table();
merged_nodule_info = [];
%% main process
for idx = 1:size(pid_list,2)
    tic
    id = strsplit(pid_list{idx},'.');
    id = id{1};
    ids = strsplit(id,'_');
    pid = ids{1}; nid = ids{2};
    nids = strsplit(nid,'-'); % to cover multiple segmentations for a nodule
    
    fprintf('%d %s %s\n', idx, pid, nid);
    %if str2num(nid) > 1, continue; end
    %if numel(strfind(pid, '_'))>0 || numel(strfind(pid, 'Ell'))>0, continue; end
    %if numel(strfind(pid, '_'))>0 || numel(strfind(pid, 'Ell'))>0 || ~strcmp(nid,'1'), continue; end
    %if numel(strfind(pid, '0043'))==0, continue; end
    

    %% spiculation quantification
    try
        spiculation_quantification
    catch
        continue
    end
    
    %% save spiculation quantification results in voxel data
    try
        voxelize_meshes
    catch
        continue
    end

    %% ignore figure generation
    continue
    
    %%
    figure('pos',[10 10 1200 600]);
    set(gcf, 'renderer', 'OpenGL');
    subplot(235), 
    bar(hx,hn/sum(hn)); hold on; 
    if numel(hx1) > 0, bar(hx1,hn1/sum(hn),'r'); alpha(0.5); end
    title('Normalized Height Distribution');
    xlim([-0.1,1.1])
    
    h1= subplot(233);
    for pki = 1:num_spikes
        %if numel(spikes(pki).baseline_curve) == 0, continue; end
        plot3(spikes(pki).baseline_curve(:,1),spikes(pki).baseline_curve(:,2),spikes(pki).baseline_curve(:,3),'r' ); hold on;
        %plot3(spikes(pki).fwhm_baseline_curve(:,1),spikes(pki).fwhm_baseline_curve(:,2),spikes(pki).fwhm_baseline_curve(:,3),'y' ); hold on;
        %plot3(spikes(pki).line_spike_base(:,1),spikes(pki).line_spike_base(:,2),spikes(pki).line_spike_base(:,3),'k'); hold on;
        plot3(spikes(pki).l_center(:,1),spikes(pki).l_center(:,2),spikes(pki).l_center(:,3),'k'); hold on;
        text(l_loc(pki,1),l_loc(pki,2),l_loc(pki,3),num2str(pki));
    end


%     patch(s, 'FaceVertexCData',t, 'FaceAlpha',0.6,'FaceColor','interp', 'EdgeColor', 'None');view(-37.5,30);
%     caxis([0 4])
%     axis equal
    hold on;
    plot3(l_loc(spiculation,1),l_loc(spiculation,2),l_loc(spiculation,3),'r*')
    plot3(l_loc(lobulation,1),l_loc(lobulation,2),l_loc(lobulation,3),'m*')
    plot3(l_loc(attached,1),l_loc(attached,2),l_loc(attached,3),'c*')

%     h2 = subplot(236);
%     copyobj(get(h1,'child'),h2);
%     %3 - default az = -37.5, el = 30.
%     view(-37.5+180,30);
%     caxis([0 4])
%     axis equal

    %
    %l_ArD=nd;
    h=subplot(232);
    [n,x]=hist(nd, 32);
    bar(x,n./sum(n),0.5); hold on;
    alpha(0.6)
    if numel([spikes.apex]) > 0
        [~,mi] = min(pdist2(x',nd([spikes.apex])));
        if numel(mi) > 0
            x1 = nd(unique(cell2mat({spikes(spiculation).vertices}')));
            [n1]=hist(x1(x1<mean(x1)), x);
            bar(x,n1./sum(n),0.5,'r','EdgeColor','none')
            alpha(0.6)
            plot(nd([spikes(spiculation).apex]),n1(mi(spiculation))./sum(n),'r*')
            %[n1]=hist(l_ArD(spikes(spiculation)), x);
            %bar(x,n1./sum(n),0.5,'r','EdgeColor','r')

            x1 = nd(unique(cell2mat({spikes(lobulation).vertices}')));
            [n1]=hist(x1(x1<mean(x1)), x);
            bar(x,n1./sum(n),0.5,'m','EdgeColor','none')
            alpha(0.6)
            plot(nd([spikes(lobulation).apex]),n1(mi(lobulation))./sum(n),'m*')
            %[n1]=hist(l_ArD(spikes(lobulation)), x);
            %bar(x,n1./sum(n),0.5,'m','EdgeColor','m')

            x1 = nd(unique(cell2mat({spikes(~(spiculation|lobulation)).vertices}')));
            [n1]=hist(x1(x1<mean(x1)), x);
            bar(x,n1./sum(n),0.5,'y','EdgeColor','none')
            alpha(0.6)
            plot(nd([spikes(~(lobulation|spiculation)).apex]),n1(mi(~(lobulation|spiculation)))./sum(n),'c*')
            %[n1]=hist(l_ArD(spikes(~(lobulation|spiculation))), x);
            %
        end
    end
    
    cmap = colormap('lines');
    xl = get(h, 'xlim');
    x1 = xl(1):0.01:xl(2);
%    y1 = gmm.pdf(x1');
    %y1 = normpdf(x1',model.muhat,model.sigmahat);
    %y1 = normpdf(x1',0,1);
    %w1 = 1/max(y1)*max(n./sum(n));
    %plot(x1,w1*y1,'Color',cmap(1,:),'LineWidth',1)
%     [~,gidxs] = sort(gmm.mu,'descend');
%     for gi = gidxs'
%         y1 = gmm.PComponents(gi)*normpdf(x1,gmm.mu(gi),sqrt(gmm.Sigma(gi)));
%         plot(x1,w1*y1,'Color',cmap(gi+1,:),'LineWidth',1)
%     end
%    y1 = normpdf(x1,muhat,sigmahat);
%    w1 = 1/max(y1)*max(n./sum(n));
%    plot(x1,w1*y1,'r--')
    xlim(xl)
    title('Area Distortion Histogram')
    
    
    subplot(231),
    %colors = (diff_area>-0.5).*sharpness;
    title('Original Mesh Model')
    patch(s, 'FaceVertexCData',colors, 'VertexNormals', normals, 'FaceAlpha',1,'FaceColor','interp', 'EdgeColor', 'interp'); hold on;
    hold on;
%     plot3(l_loc(:,1),l_loc(:,2),l_loc(:,3),'y*')
    plot3(l_loc(spiculation,1),l_loc(spiculation,2),l_loc(spiculation,3),'r*')
    plot3(l_loc(lobulation,1),l_loc(lobulation,2),l_loc(lobulation,3),'m*')
    plot3(l_loc(attached,1),l_loc(attached,2),l_loc(attached,3),'c*')

    view(-37.5,30);
    caxis([-max(abs(l_ArD)) max(abs(l_ArD))])
    axis equal fill
    
   
%    subplot(232),
%     title('Area Distortion Map')
%     %colors = (diff_area>-0.5).*sharpness;
%     %colors = l_ArD;
%     %patch(s, 'FaceAlpha',1,'FaceColor','interp', 'EdgeColor', 'interp');hold on
%     %patch(s,'FaceColor','red','FaceAlpha',0.5,'EdgeColor','none'); hold on
%     patch(s, 'FaceVertexCData',l_ArD, 'VertexNormals', normals, 'FaceAlpha',0.9,'EdgeAlpha',0.9,'FaceColor','interp', 'EdgeColor', 'interp');
%     hold on;
%     colorbar('eastoutside')
%     plot3(l_loc(:,1),l_loc(:,2),l_loc(:,3),'y*')
%     plot3(l_loc(spiculation,1),l_loc(spiculation,2),l_loc(spiculation,3),'r*')
%     plot3(l_loc(lobulation,1),l_loc(lobulation,2),l_loc(lobulation,3),'m*')
% 
%     view(-37.5,30);
%     caxis([-max(abs(l_ArD)) max(abs(l_ArD))])
%     axis equal fill
%     %camlight;  camlight(-80,-10); lighting phong;

    colormap(jet)

    h1 = subplot(233);
    patch(s, 'FaceVertexCData',nd,'FaceAlpha',0.6,'FaceColor','interp', 'EdgeColor', 'none');view(-37.5,30);
    %patch(s, 'FaceVertexCData',d,'FaceAlpha',0.8,'FaceColor','interp', 'EdgeColor', 'interp');view(-37.5,30);
    hold on;
%     plot3(l_loc(:,1),l_loc(:,2),l_loc(:,3),'y*')
    plot3(l_loc(spiculation,1),l_loc(spiculation,2),l_loc(spiculation,3),'r*')
    plot3(l_loc(lobulation,1),l_loc(lobulation,2),l_loc(lobulation,3),'m*')
    plot3(l_loc(attached,1),l_loc(attached,2),l_loc(attached,3),'c*')
    %patch('Vertices', nv, 'Faces', s.faces, 'FaceVertexCData',Cgaussian, 'FaceColor','interp', 'EdgeColor', 'interp');view(-37.5,30);
    %caxis([0 5])
    caxis([-1 1])
    view(-37.5,30);
    colorbar('eastoutside', 'YAxisLocation', 'right');
    axis equal fill
    title('spike detection (front)')
    %camlight;  camlight(-80,-10); lighting phong;


    subplot(234),
    %colors = (diff_area>-0.5).*sharpness;
    %colors = l_ArD;
    l_loc1 = s1.vertices([spikes.apex],:);
    %patch(s1, 'FaceAlpha',1,'FaceColor','green', 'EdgeColor', 'interp'); hold on
    %patch(s1,'FaceColor','red','FaceAlpha',0.5,'EdgeColor','none'); hold on
    title('Spherical Mapping')
    patch(s1, 'FaceVertexCData',colors, 'VertexNormals', normals, 'FaceAlpha',1,'FaceColor','interp', 'EdgeColor', 'interp'); hold on;
%     plot3(l_loc1(:,1),l_loc1(:,2),l_loc1(:,3),'y*')
    plot3(l_loc1(spiculation,1),l_loc1(spiculation,2),l_loc1(spiculation,3),'r*')
    plot3(l_loc1(lobulation,1),l_loc1(lobulation,2),l_loc1(lobulation,3),'m*')
    plot3(l_loc1(attached,1),l_loc1(attached,2),l_loc1(attached,3),'c*')
    view(-37.5,30);
    axis equal fill


    h2 = subplot(236);
    for pki = 1:num_spikes
        %if numel(spikes(pki).baseline_curve) == 0, continue; end
        plot3(spikes(pki).baseline_curve(:,1),spikes(pki).baseline_curve(:,2),spikes(pki).baseline_curve(:,3),'r' ); hold on;
        %plot3(spikes(pki).fwhm_baseline_curve(:,1),spikes(pki).fwhm_baseline_curve(:,2),spikes(pki).fwhm_baseline_curve(:,3),'y' ); hold on;
        %plot3(spikes(pki).line_spike_base(:,1),spikes(pki).line_spike_base(:,2),spikes(pki).line_spike_base(:,3),'k'); hold on;
        plot3(spikes(pki).l_center(:,1),spikes(pki).l_center(:,2),spikes(pki).l_center(:,3),'k'); hold on;
        text(l_loc(pki,1),l_loc(pki,2),l_loc(pki,3),num2str(pki));
    end
%     patch(s1, 'FaceVertexCData',l_ArD,'FaceAlpha',0.8,'FaceColor','interp', 'EdgeColor', 'interp');view(-37.5,30);
%     hold on; plot3(l_loc(:,1),l_loc(:,2),l_loc(:,3),'r*')
    %patch('Vertices', nv, 'Faces', s.faces, 'FaceVertexCData',Cgaussian, 'FaceColor','interp', 'EdgeColor', 'interp');view(-37.5,30);
    patch(s, 'FaceVertexCData',nd,'FaceAlpha',0.6,'FaceColor','interp', 'EdgeColor', 'none');view(-37.5,30);
    %patch(s, 'FaceVertexCData',d,'FaceAlpha',0.8,'FaceColor','interp', 'EdgeColor', 'interp');view(-37.5,30);
    hold on;
%     plot3(l_loc(:,1),l_loc(:,2),l_loc(:,3),'y*')
    plot3(l_loc(spiculation,1),l_loc(spiculation,2),l_loc(spiculation,3),'r*')
    plot3(l_loc(lobulation,1),l_loc(lobulation,2),l_loc(lobulation,3),'m*')
    plot3(l_loc(attached,1),l_loc(attached,2),l_loc(attached,3),'c*')
    %patch('Vertices', nv, 'Faces', s.faces, 'FaceVertexCData',Cgaussian, 'FaceColor','interp', 'EdgeColor', 'interp');view(-37.5,30);
    %caxis([0 5])
    caxis([-1 1])
    view(-37.5+180,30);
    colorbar('eastoutside', 'YAxisLocation', 'right');
    axis equal fill
    title('spike detection (back)')
    
    
    
    title_text = [id,' S', num2str(n_info.Spiculation),' M',num2str(n_info.Malignancy),' Volume:',num2str(volume/1000,'%0.2f'),'cc, D:',num2str(esr*2,'%0.2f'), ...
        'mm, Roundness:',num2str(roundness,'%0.2f'),', #Spikes:',num2str(num_spikes),', #Spiculation:',num2str(sum(spiculation)),', S1:',num2str(spic1,'%0.2f'),', S2:',num2str(spic2,'%0.2f'), ...
        ', S3:',num2str(spic3,'%0.2f'),', S4:',num2str(spic4,'%0.2f')];
    %titile_text = [id,' Volume:',num2str(volume/1000,2),'cc, D:',num2str(esr*2),'mm, #Spikes:',num2str(num_spikes),', Roundness:',num2str(roundness,2)];
    h = sgtitle(title_text); set(h, 'FontSize',12,'Interpreter','none');
    

    %%
    %writetable(f, ['output/' pid '/' pid '_' nid '_' type num2str(sigma,'%.1f') '-' num2str(gap) '-' num2str(is_boundary_only) num2str(is_contour) num2str(weight_mode) '.csv']);
    
    savefig([output_experiment_path '/spiculation/' pid '_' nid '_sphere_param.fig'])

    set(gcf, 'PaperUnits', 'inches');
    x_width=12; y_width=6;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


    print([output_experiment_path '/spiculation/' pid '_' nid '_sphere_param.png'],'-dpng','-r300')
    close all
    toc
end
writetable(all_spikes, [output_experiment_path '/spiculation/spikes.csv']);
writetable(features, [output_experiment_path '/spiculation/features.csv']);
