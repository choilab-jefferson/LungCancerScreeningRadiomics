function [f, spikes_table]  = spiculation_pipeline(pid, nid, filename, n_info, ...
    label_file, attached_label_filename, obj_filename, sph_map_filename, ...
    ard_voxel_filename, ard_surface_filename, spike_lable_filename, spike_surface_filename, ...
    ard_filename, attached_indices_filename, spike_apex_filename, spiculation_fcsv_filename, feature_filename, ...
    figure_fig_filename, figure_png_filename)
    %%
    try
        [s.vertices,s.faces] = readOBJ(obj_filename);
        try
            [s1.vertices,s1.faces,~,normals] = read_obj_colors(sph_map_filename);
        catch
            [s1.vertices,s1.faces] = readOBJ(sph_map_filename);
            if size(s1.vertices,2) == 6
                colors = s1.vertices(:,4:6);
                s1.vertices = s1.vertices(:,1:3);
            end
            [normals,normals_f] = patchnormals(s1);
        end
        assert(numel(s.vertices)==numel(s1.vertices));
    catch exception
        if ~exist(obj_filename, 'file')
            display('no obj file')
        end
        if ~exist(sph_map_filename, 'file') || numel(s.vertices)~=numel(s1.vertices)
            display('no spherical mapping!')
            s1=[];
            normals=[];
            %[s1,~,normals] = sphereical_mapping(s);
        end
        throw(exception)
    end
    if ~exist('colors', 'var')
        [normal_s,normalf_s] = patchnormals(s);
        colors = (1-normal_s)/2; % recaculate vetex color
        colors = 1-colors; % flip color space to show bright color
    end

    %% attachment vertices
    [o_attached_img_3d, meta] = fn_nrrdread(attached_label_filename);
    attachment = find_attached_vertices(o_attached_img_3d, meta, s);
    attachment_indices = find(attachment);


    %% spiculation quantification
    [spikes, volume, area, nd, l_ArD, l_AnD] = spiculation_quantification(s, s1, attachment_indices);

    %%
    %figure, hist([spikes.height])

    %% extract shape features from nodule mesh model
    esph_factor = (3/(4*pi))^(1/3);
    esr = esph_factor*volume^(1/3);
    roundness = min((4*pi*esr^2)/area,1);

    [num_spikes,l_height,l_width,l_angle,l_nheight,l_nwidth,l_rimchange, ...
        spic_ap,spic_bp,spic1p,spic2p,spic3p,spic4p,hx,hn] = extract_spiculation_features(spikes, esr);
    nhn = hn / sum(hn);


    %% convert spikes structure to table
    if num_spikes > 0
        if num_spikes > 1
            spikes_table = struct2table(spikes);
        else
            spikes(2) = spikes(1);
            spikes_table = struct2table(spikes);
            spikes_table = spikes_table(1,:);
            spikes(2) = [];
        end
        spikes_table.Filename = repmat({filename},num_spikes,1);
        spikes_table.PID = repmat({pid},num_spikes,1);
        spikes_table.NID = repmat({nid},num_spikes,1);
        spikes_table.D = repmat(esr*2,num_spikes,1);
        spikes_table.volume = repmat(volume/1000,num_spikes,1);
        spikes_table.nd = nd(spikes_table.apex);
        spikes_table.rimchange = l_rimchange';
        spikes_table.nrimchange = l_rimchange'./(l_height'+1);
        if numel(n_info) > 0
            spikes_table.spic = repmat(n_info.Spiculation,num_spikes,1);
            spikes_table.lobul = repmat(n_info.Lobulation,num_spikes,1);
            spikes_table.malig = repmat(n_info.Malignancy,num_spikes,1);
            spikes_table = spikes_table(:,{'Filename','PID','NID','spic','lobul','malig','volume','D','nd','nrimchange','apex','height','width','angle'});
        else
            spikes_table = spikes_table(:,{'Filename','PID','NID','volume','D','nd','nrimchange','apex','height','width','angle'});
        end

        disp(spikes_table)
    end


    %% remove small spikes
    th_noise = 3;
    selected_spikes = l_height > th_noise & l_width > th_noise/2;
    spikes = spikes(selected_spikes);
    num_spikes = numel(spikes);
    fprintf('number of spikes %d\n', num_spikes)

    %%
    l_loc = s.vertices([spikes.apex],:); % list of apex coordinates


    %% spike classification, spiculation or lobulation
    spiculation = [spikes.type]==0;
    lobulation = [spikes.type]==1;
    attached = [spikes.type]==3;


    %% recalculate features
    [num_spikes,l_height,l_width,l_angle,l_nheight,l_nwidth,l_rimchange, ...
        spic_a,spic_b,spic1,spic2,spic3,spic4,hx1,hn1] = extract_spiculation_features(spikes, esr);
    nhn1 = hn1 / sum(hn1);


    %% save important parameters
    if exist('ard_filename', 'var')
        dlmwrite(ard_filename, nd);
    end
    if exist('attached_indices_filename', 'var')
        dlmwrite(attached_indices_filename, attachment_indices);
    end
    if exist('spike_apex_filename', 'var')
        dlmwrite(spike_apex_filename,[spikes.apex]');
    end

    %% save apex point of spikes and area distortion for slicer
    if exist('spiculation_fcsv_filename', 'var')
        fid = fopen(spiculation_fcsv_filename,'w');
        fprintf(fid, '# Markups fiducial file version = 4.8\n');
        fprintf(fid, '# CoordinateSystem = 0\n');
        fprintf(fid, '# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n');
        for pki = 1:num_spikes
            if spikes(pki).type == 0
                fprintf(fid, 'vtkMRMLMarkupsFiducialNode_%d,%f,%f,%f,0,0,0,0,1,1,1,S%d,"AD %0.2f",vtkMRMLScalarVolumeNode1\n',pki,l_loc(pki,1),l_loc(pki,2),l_loc(pki,3),pki,nd(spikes(pki).apex));
            elseif spikes(pki).type == 1
                fprintf(fid, 'vtkMRMLMarkupsFiducialNode_%d,%f,%f,%f,0,0,0,0,1,0,2,L%d,"AD %0.2f",vtkMRMLScalarVolumeNode1\n',pki,l_loc(pki,1),l_loc(pki,2),l_loc(pki,3),pki,nd(spikes(pki).apex));
            elseif spikes(pki).type == 2
                fprintf(fid, 'vtkMRMLMarkupsFiducialNode_%d,%f,%f,%f,0,0,0,0,0,0,2,C%d,"AD %0.2f",vtkMRMLScalarVolumeNode1\n',pki,l_loc(pki,1),l_loc(pki,2),l_loc(pki,3),pki,nd(spikes(pki).apex));
            elseif spikes(pki).type == 3
                fprintf(fid, 'vtkMRMLMarkupsFiducialNode_%d,%f,%f,%f,0,0,0,0,0,0,2,A%d,"AD %0.2f",vtkMRMLScalarVolumeNode1\n',pki,l_loc(pki,1),l_loc(pki,2),l_loc(pki,3),pki,nd(spikes(pki).apex));
            end
        end
        fclose(fid);
    end


    %% raw features from curvature maps
    [Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2]=patchcurvature(s);

    l_GC = Cgaussian(:);
    l_MC = Cmean(:); 

    %% spike detection by mean curvature
    % spikes_old = spiculation_detection_mean_curvature(Cmean, Cgaussian, s, s1, nd); % baseline detection - area distortion still required


    %% feature table
    hnames = strsplit(sprintf('h%d ',1:10));
    nhnames = strsplit(sprintf('nh%d ',1:10));

    f = table({filename}, {pid},{nid},num_spikes,sum(spiculation),sum(lobulation), sum(attached), mean(attachment), ...
            volume,area,esr*2,roundness,spic_a, spic_b, spic_ap, spic_bp, ...
            spic1,spic2,spic3,spic4,spic1p,spic2p,spic3p,spic4p,...
            min(l_GC),max(l_GC),median(l_GC),mean(l_GC),var(l_GC),skewness(l_GC),kurtosis(l_GC), ...
            min(l_MC), max(l_MC),median(l_MC),mean(l_MC),var(l_MC),skewness(l_MC),kurtosis(l_MC), ...
            min(nd),max(nd),median(nd),mean(nd),var(nd),skewness(nd),kurtosis(nd), ...
            min(l_ArD),max(l_ArD),median(l_ArD),mean(l_ArD),var(l_ArD),skewness(l_ArD),kurtosis(l_ArD), ...
            min(l_AnD), max(l_AnD),median(l_AnD),mean(l_AnD),var(l_AnD),skewness(l_AnD),kurtosis(l_AnD), ...
            min(l_height), max(l_height),median(l_height),mean(l_height),var(l_height),skewness(l_height),kurtosis(l_height), ...
            min(l_width), max(l_width),median(l_width),mean(l_width),var(l_width),skewness(l_width),kurtosis(l_width), ...
            min(l_angle), max(l_angle),median(l_angle),mean(l_angle),var(l_angle),skewness(l_angle),kurtosis(l_angle), ...
            min(l_nheight), max(l_nheight),median(l_nheight),mean(l_nheight),var(l_nheight),skewness(l_nheight),kurtosis(l_nheight), ...
            min(l_nwidth), max(l_nwidth),median(l_nwidth),mean(l_nwidth),var(l_nwidth),skewness(l_nwidth),kurtosis(l_nwidth), ...
            hn(1),hn(2),hn(3),hn(4),hn(5),hn(6),hn(7),hn(8),hn(9),hn(10), ...
            nhn(1),nhn(2),nhn(3),nhn(4),nhn(5),nhn(6),nhn(7),nhn(8),nhn(9),nhn(10), ...
            'VariableNames',{'Filename','PID','NID','num_spikes','num_spic','num_lobul','num_attached','attachment', ...
            'volume','area','eqv_sph_D','roundness','spic_a','spic_b','spic_ap','spic_bp', ...
            'spic1','spic2','spic3','spic4','spic1p','spic2p','spic3p','spic4p', ...
            'min_GC','max_GC','med_GC','mean_GC','var_GC','skew_GC','kurt_GC', ...
            'min_MC','max_MC','med_MC','mean_MC','var_MC','skew_MC','kurt_MC', ...
            'min_NArD','max_NArD','med_NArD','mean_NArD','var_NArD','skew_NArD','kurt_NArD', ...
            'min_ArD','max_ArD','med_ArD','mean_ArD','var_ArD','skew_ArD','kurt_ArD', ...
            'min_AnD','max_AnD','med_AnD','mean_AnD','var_AnD','skew_AnD','kurt_AnD', ...
            'min_h','max_h','med_h','mean_h','var_h','skew_h','kurt_h', ...
            'min_w','max_w','med_w','mean_w','var_w','skew_w','kurt_w', ...
            'min_a','max_a','med_a','mean_a','var_a','skew_a','kurt_a', ...
            'min_nh','max_nh','med_nh','mean_nh','var_nh','skew_nh','kurt_nh', ...
            'min_nw','max_nw','med_nw','mean_nw','var_nw','skew_nw','kurt_nw' ...
            hnames{1:10},nhnames{1:10}});
    if exist('feature_filename', 'var')
        writetable(f, feature_filename);
    end
 


    %% save spiculation quantification results in voxel data
    try
        [o_seg_img_3d, meta] = fn_nrrdread(label_file);
        surface_voxel = bwperim(o_seg_img_3d, 6);
        sigma = esr/2/3;
        
        [ard_voxel, spike_label] = voxelize_meshes(o_seg_img_3d, meta, s, nd, spikes, sigma);
        
        fn_nrrdwrite(ard_voxel_filename, ard_voxel, meta);
        fn_nrrdwrite(ard_surface_filename, ard_voxel.*surface_voxel, meta);
        fn_nrrdwrite(spike_lable_filename, spike_label, meta);
        fn_nrrdwrite(spike_surface_filename, spike_label.*surface_voxel, meta);
    catch exception
        throw(exception)
    end


    %% Figure
    if ~exist('figure_fig_filename', 'var') and ~exist('figure_fig_filename', 'var')
        return
    end

    %%
    figure('pos',[10 10 1200 600]);
    %set(gcf, 'renderer', 'OpenGL');
    subplot(235), 
    bar(hx,hn/sum(hn)); hold on; 
    if numel(hx1) > 0, bar(hx1,hn1/sum(hn),'r'); alpha(0.5); end
    title('Normalized Height Distribution');
    xlim([-0.1,1.1])

    h1= subplot(233);
    for pki = 1:num_spikes
        %if numel(spikes(pki).baseline_curve) == 0, continue; end
        plot3(spikes(pki).baseline_curve(:,1),spikes(pki).baseline_curve(:,2),spikes(pki).baseline_curve(:,3),'r' ); hold on;
        plot3(spikes(pki).fwhm_baseline_curve(:,1),spikes(pki).fwhm_baseline_curve(:,2),spikes(pki).fwhm_baseline_curve(:,3),'y' ); hold on;
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
            plot(nd([spikes(~(lobulation|spiculation)).apex]),n1(mi(~(lobulation|spiculation)))./sum(n),'y*')
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
    view(-37.5,30);
    axis equal fill


    h2 = subplot(236);
    for pki = 1:num_spikes
        %if numel(spikes(pki).baseline_curve) == 0, continue; end
        plot3(spikes(pki).baseline_curve(:,1),spikes(pki).baseline_curve(:,2),spikes(pki).baseline_curve(:,3),'r' ); hold on;
        plot3(spikes(pki).fwhm_baseline_curve(:,1),spikes(pki).fwhm_baseline_curve(:,2),spikes(pki).fwhm_baseline_curve(:,3),'y' ); hold on;
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
    %patch('Vertices', nv, 'Faces', s.faces, 'FaceVertexCData',Cgaussian, 'FaceColor','interp', 'EdgeColor', 'interp');view(-37.5,30);
    %caxis([0 5])
    caxis([-1 1])
    view(-37.5+180,30);
    colorbar('eastoutside', 'YAxisLocation', 'right');
    axis equal fill
    title('spike detection (back)')



    %title_text = [id,' S', num2str(n_info.Spiculation),' M',num2str(n_info.Malignancy),' Volume:',num2str(volume/1000,'%0.2f'),'cc, D:',num2str(esr*2,'%0.2f'), ...
    %    'mm, Roundness:',num2str(roundness,'%0.2f'),', #Spikes:',num2str(num_spikes),', S1:',num2str(spic1,'%0.2f'),', S1:',num2str(spic2,'%0.2f'), ...
    %    ', S3:',num2str(spic3,'%0.2f'),', S4:',num2str(spic4,'%0.2f')];
    title_text = [pid, '_', nid,' Volume:',num2str(volume/1000,2),'cc, D:',num2str(esr*2),'mm, #Spikes:',num2str(num_spikes),', #Spiculation:',num2str(sum(spiculation)),', Roundness:',num2str(roundness,2)];
    h = sgtitle(title_text); set(h, 'FontSize',12,'Interpreter','none');


    %%
    %writetable(f, ['output/' pid '/' pid '_' nid '_' type num2str(sigma,'%.1f') '-' num2str(gap) '-' num2str(is_boundary_only) num2str(is_contour) num2str(weight_mode) '.csv']);

    if exist('figure_fig_filename', 'var')
        savefig(figure_fig_filename)
    end

    set(gcf, 'PaperUnits', 'inches');
    x_width=12; y_width=6;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

    if exist('figure_png_filename', 'var')
        print(figure_png_filename,'-dpng','-r300')
    end
end
