function [rt_N, vt_N] = observatoryPosition(t)
%OBSERVATORYPOSITION Summary of this function goes here
%   Assume prime meridian in inertial X direction at t=0

latitude = 28.467*pi/180; % rad
longitude = -80.4767*pi/180; % rad

T_earth = 23.9345; %[h]
R_earth = 6378; % [km]
omega_earth = 2*pi/(T_earth*60^2); %[rad/s]

phi = wrapTo180(latitude + t*omega_earth);

rt_N = Euler3(phi)*Euler2(longitude)*[R_earth; 0; 0];

vt_N = Euler3(phi)*[0; R_earth*omega_earth; 0];
end

