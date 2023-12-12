function [sigRN, omegaRN_R, P] = missionTracking(t)
    
    % epoch0 = 05/23/1999 00:16:12.24
    rs_N = [1;0;0]; % Assume sun position is constant in inertial X direction    
    vs_N = [0;0;0];
    rs = norm(rs_N);
    
    [rt_N, vt_N] = observatoryPosition(t);
    [rc_N, vc_N] = orbitPropogator(t);
    
    ls_N = rs_N - rc_N;
    ls = norm(ls_N);
    lt_N = rt_N - rc_N;
    lt = norm(lt_N);
    
    d_ls_N = vs_N - vc_N;
    d_lt_N = vt_N - vc_N;
    
    % want z axis alligned with lt and y axis perpendicular to ls
    
    z_N = lt_N/lt;
    y_N = cross(z_N,ls_N)/norm(cross(z_N,ls_N));
    x_N = cross(y_N, z_N);
    
    RN = [x_N'; y_N'; z_N'];
    sigRN = C2MRP(RN);
    
    omegaX = -dot(d_lt_N, y_N)/lt;
    omegaY = dot(d_lt_N, x_N)/lt;
    omegaZ = (dot(omegaX*ls_N,z_N) + dot(y_N, d_ls_N)) / dot(ls_N, x_N);
    
    omegaRN_N = -[omegaX; omegaY; omegaZ];
    omegaRN_R = RN*omegaRN_N;
    
    % booleon for if the spacecraft is in eclipse
    R_earth = 6378; % [km]
    eclipse = (dot(rc_N, rs_N) < 0) & (norm(cross(rc_N, rs_N/rs)) < R_earth);
    
    P_eclipse = -680; %[W] wuthdraw during eclipse
    P_sun = 1000; %[W]
    if eclipse == 1
        P = P_eclipse;
    else
        P = P_sun;
    end
    
end

