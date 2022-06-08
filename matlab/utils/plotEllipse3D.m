function [rotated_ellipse] = plotEllipse3D(center,normal,ellipse_t)
    a = ellipse_t.a;
    b = ellipse_t.b;
    X0 = ellipse_t.X0;
    Y0 = ellipse_t.Y0;
    cos_phi = cos( ellipse_t.phi );
    sin_phi = sin( ellipse_t.phi );

    % rotation matrix to rotate the axes with respect to an angle phi
    R = [ cos_phi sin_phi 0; -sin_phi cos_phi 0; 0 0 1 ];
    
    % the axes
    ver_line        = [ [X0 X0]; Y0+b*[-1 1]; 0 0];
    horz_line       = [ X0+a*[-1 1]; [Y0 Y0]; 0 0];
    new_ver_line    = R*ver_line;
    new_ver_line    = RodriguesRotation(new_ver_line',[0, 0, 1],normal')'+center';
    new_horz_line   = R*horz_line;
    new_horz_line   = RodriguesRotation(new_horz_line',[0, 0, 1],normal')'+center';
    
    % the ellipse
    theta_r         = linspace(0,2*pi);
    ellipse_x_r     = X0 + a*cos( theta_r );
    ellipse_y_r     = Y0 + b*sin( theta_r );
    rotated_ellipse = R * [ellipse_x_r;ellipse_y_r; zeros(size(ellipse_y_r))];
    rotated_ellipse = RodriguesRotation(rotated_ellipse', [0, 0, 1], normal')'+center';

    
    % draw
    %plot3( new_ver_line(1,:),new_ver_line(2,:),new_ver_line(3,:),'r' ); hold on
    %plot3( new_horz_line(1,:),new_horz_line(2,:),new_horz_line(3,:),'r' );
    plot3( rotated_ellipse(1,:),rotated_ellipse(2,:),rotated_ellipse(3,:),'r' );
end