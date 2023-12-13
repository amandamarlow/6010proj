function [rt_N, vt_N] = observatoryPosition(t)
%OBSERVATORYPOSITION Summary of this function goes here
%   Assume prime meridian in inertial X direction at t=0

% latitude = 28.467*pi/180; % rad
% longitude = -80.4767*pi/180; % rad

T_earth = 23.9345; %[h]
% R_earth = 6373.3; % [km] % radius at latitude
omega_earth = 2*pi/(T_earth*60^2); %[rad/s]

% phi = wrapTo180(longitude + t*omega_earth);

r0_N = [-5386.525298; 1572.798194; 3021.709804];
% R_earth = norm(r0_N);
theta = t*omega_earth;
rt_N = Euler3(theta)*r0_N;

% rt_N = (Euler3(phi)*Euler2(latitude))'*[R_earth; 0; 0];

vt_N = (Euler3(theta))'*[0; norm(r0_N(1:2))*omega_earth; 0];
end

