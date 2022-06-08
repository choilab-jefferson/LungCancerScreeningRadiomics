%% main function
clc;
clear;
close all

set_environment; % call environment setting

%% spiculation toolbox
addpath(genpath('utils'))
addpath(genpath('spiculation_detection'))

%% set values
experiment_set = 'nodule-lungx';
esph_factor = (3/(4*pi))^(1/3);
subset = '';%'_smooth';

%% directory paths
experiment_path = [data_path '/' experiment_set];
output_experiment_path=[output_path '/' experiment_set];
if ~isdir(output_experiment_path); mkdir(output_experiment_path); end
if ~isdir([output_experiment_path '/parameters']); mkdir([output_experiment_path '/parameters']); end

obj_path = [experiment_path '/objs' subset];
sph_map_path = [experiment_path '/spherical_obj' subset];
dir_obj = dir([obj_path '/*.obj']);
pid_list = {dir_obj.name};

features = table();
all_peaks = table();
merged_nodule_info = [];
%% main process
for idx = 1:size(pid_list,2)
    tic % tic starts a stopwatch timer

    %% extract PID and NID from file name
    id = strsplit(pid_list{idx},'.');
    id = id{1};
    ids = strsplit(id,'_');
    pid = ids{1}; nid = ids{2};
    nids = strsplit(nid,'-');

    fprintf('%d %s %s\n', idx, pid, nid);

    %% ignore some cases
    if numel(strfind(pid, '_'))>0 || numel(strfind(pid, 'Ell'))>0, continue; end
    %if numel(strfind(pid, '_'))>0 || numel(strfind(pid, 'Ell'))>0 || ~strcmp(nid,'1'), continue; end    
    %if numel(strfind(pid, '0289'))==0, continue; end
    
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
    for pki = 1:num_peaks
        %if numel(peaks(pki).baseline_curve) == 0, continue; end
        plot3(peaks(pki).baseline_curve(:,1),peaks(pki).baseline_curve(:,2),peaks(pki).baseline_curve(:,3),'r' ); hold on;
        plot3(peaks(pki).fwhm_baseline_curve(:,1),peaks(pki).fwhm_baseline_curve(:,2),peaks(pki).fwhm_baseline_curve(:,3),'y' ); hold on;
        %plot3(peaks(pki).line_peak_base(:,1),peaks(pki).line_peak_base(:,2),peaks(pki).line_peak_base(:,3),'k'); hold on;
        plot3(peaks(pki).l_center(:,1),peaks(pki).l_center(:,2),peaks(pki).l_center(:,3),'k'); hold on;
        text(l_loc(pki,1),l_loc(pki,2),l_loc(pki,3),num2str(pki));
    end


%     patch(s, 'FaceVertexCData',t, 'FaceAlpha',0.6,'FaceColor','interp', 'EdgeColor', 'None');view(-37.5,30);
%     caxis([0 4])
%     axis equal
    hold on;
    plot3(l_loc(spiculation,1),l_loc(spiculation,2),l_loc(spiculation,3),'r*')
    plot3(l_loc(lobulation,1),l_loc(lobulation,2),l_loc(lobulation,3),'m*')

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
    if numel([peaks.apex]) > 0
        [~,mi] = min(pdist2(x',nd([peaks.apex])));
        if numel(mi) > 0
            x1 = nd(unique(cell2mat({peaks(spiculation).vertices}')));
            [n1]=hist(x1(x1<mean(x1)), x);
            bar(x,n1./sum(n),0.5,'r','EdgeColor','none')
            alpha(0.6)
            plot(nd([peaks(spiculation).apex]),n1(mi(spiculation))./sum(n),'r*')
            %[n1]=hist(l_ArD(peaks(spiculation)), x);
            %bar(x,n1./sum(n),0.5,'r','EdgeColor','r')

            x1 = nd(unique(cell2mat({peaks(lobulation).vertices}')));
            [n1]=hist(x1(x1<mean(x1)), x);
            bar(x,n1./sum(n),0.5,'m','EdgeColor','none')
            alpha(0.6)
            plot(nd([peaks(lobulation).apex]),n1(mi(lobulation))./sum(n),'m*')
            %[n1]=hist(l_ArD(peaks(lobulation)), x);
            %bar(x,n1./sum(n),0.5,'m','EdgeColor','m')

            x1 = nd(unique(cell2mat({peaks(~(spiculation|lobulation)).vertices}')));
            [n1]=hist(x1(x1<mean(x1)), x);
            bar(x,n1./sum(n),0.5,'y','EdgeColor','none')
            alpha(0.6)
            plot(nd([peaks(~(lobulation|spiculation)).apex]),n1(mi(~(lobulation|spiculation)))./sum(n),'y*')
            %[n1]=hist(l_ArD(peaks(~(lobulation|spiculation))), x);
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
    %patch('Vertices', nv, 'Faces', s.faces, 'FaceVertexCData',Cgaussian, 'FaceColor','interp', 'EdgeColor', 'interp');view(-37.5,30);
    %caxis([0 5])
    caxis([-1 1])
    view(-37.5,30);
    colorbar('eastoutside', 'YAxisLocation', 'right');
    axis equal fill
    title('Peak detection (front)')
    %camlight;  camlight(-80,-10); lighting phong;


    subplot(234),
    %colors = (diff_area>-0.5).*sharpness;
    %colors = l_ArD;
    l_loc1 = s1.vertices([peaks.apex],:);
    %patch(s1, 'FaceAlpha',1,'FaceColor','green', 'EdgeColor', 'interp'); hold on
    %patch(s1,'FaceColor','red','FaceAlpha',0.5,'EdgeColor','none'); hold on
    title('Spherical Mapping')
    patch(s1, 'FaceVertexCData',colors, 'VertexNormals', normals, 'FaceAlpha',1,'FaceColor','interp', 'EdgeColor', 'interp'); hold on;
%     plot3(l_loc1(:,1),l_loc1(:,2),l_loc1(:,3),'y*')
    plot3(l_loc1(spiculation,1),l_loc1(spiculation,2),l_loc1(spiculation,3),'r*')
    plot3(l_loc1(lobulation,1),l_loc1(lobulation,2),l_loc1(lobulation,3),'m*')
    view(-37.5,30);
    axis equal fill


    h2 = subplot(236);
    for pki = 1:num_peaks
        %if numel(peaks(pki).baseline_curve) == 0, continue; end
        plot3(peaks(pki).baseline_curve(:,1),peaks(pki).baseline_curve(:,2),peaks(pki).baseline_curve(:,3),'r' ); hold on;
        plot3(peaks(pki).fwhm_baseline_curve(:,1),peaks(pki).fwhm_baseline_curve(:,2),peaks(pki).fwhm_baseline_curve(:,3),'y' ); hold on;
        %plot3(peaks(pki).line_peak_base(:,1),peaks(pki).line_peak_base(:,2),peaks(pki).line_peak_base(:,3),'k'); hold on;
        plot3(peaks(pki).l_center(:,1),peaks(pki).l_center(:,2),peaks(pki).l_center(:,3),'k'); hold on;
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
    %patch('Vertices', nv, 'Faces', s.faces, 'FaceVertexCData',Cgaussian, 'FaceColor','interp', 'EdgeColor', 'interp');view(-37.5,30);
    %caxis([0 5])
    caxis([-1 1])
    view(-37.5+180,30);
    colorbar('eastoutside', 'YAxisLocation', 'right');
    axis equal fill
    title('Peak detection (back)')



    %title_text = [id,' S', num2str(n_info.Spiculation),' M',num2str(n_info.Malignancy),' Volume:',num2str(volume/1000,'%0.2f'),'cc, D:',num2str(esr*2,'%0.2f'), ...
    %    'mm, Roundness:',num2str(roundness,'%0.2f'),', #Peaks:',num2str(num_peaks),', S1:',num2str(spic1,'%0.2f'),', S1:',num2str(spic2,'%0.2f'), ...
    %    ', S3:',num2str(spic3,'%0.2f'),', S4:',num2str(spic4,'%0.2f')];
    title_text = [pid , nid,' Volume:',num2str(volume/1000,2),'cc, D:',num2str(esr*2),'mm, #Peaks:',num2str(num_peaks),', #Spiculation:',num2str(sum(spiculation)),', Roundness:',num2str(roundness,2)];
    h = sgtitle(title_text); set(h, 'FontSize',12,'Interpreter','none');


    %%
    %writetable(f, ['output/' pid '/' pid '_' nid '_' type num2str(sigma,'%.1f') '-' num2str(gap) '-' num2str(is_boundary_only) num2str(is_contour) num2str(weight_mode) '.csv']);

    savefig([output_experiment_path '/' pid '_' nid '_sphere_param.fig'])

    set(gcf, 'PaperUnits', 'inches');
    x_width=12; y_width=6;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


    print([output_experiment_path '/' pid '_' nid '_sphere_param.png'],'-dpng','-r300')
    close all
    toc
end
writetable(all_peaks, [output_experiment_path '/peaks.csv']);
writetable(features, [output_experiment_path '/features.csv']);
