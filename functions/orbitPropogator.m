function [rc_N, vc_N] = orbitPropogator(t)
%ORBITPROPOGATOR Summary of this function goes here
%   Detailed explanation goes here

mu = 398600; %km^3/s2

% given in "Satellite Attitude Control and Power Tracking
% with Energy/Momentum Wheels"
n = 14.57788549/24/60^2*2*pi; % [rad/s]
M0 = 234.7460*pi/180; % [rad] mean anomaly
omega = 125.5766*pi/180; % [rad]
OMEGA = 132.8782*pi/180; % [rad]
i = 86.5318*pi/180; % [rad]
e = 0.00216220; 

a = (mu/n^2)^(1/3);

E0 = keplers(M0, a, e, mu);
theta0 = 2*atan(sqrt((1+e)/(1-e))*tan(E0/2));
% r0 = a*(1-e*cos(E0));
p = a*(1-e^2);
r0 = p/(1+e*cos(theta0));
% theta0 = acos((p/r0 - 1)/e);
r0_P = r0*[cos(theta0); sin(theta0); 0]; % initial position in perifocal frame
h = sqrt(mu*p);
v0_P = mu/h*[-sin(theta0); (e+cos(theta0)); 0];

Mt = M0 + n*t;
Et = keplers(Mt, a, e, mu);

[rc_P, vc_P] = fg(Et, E0, r0_P, v0_P, t, mu);

PN = Euler3(OMEGA)*Euler1(i)*Euler3(omega);
NP = PN';
rc_N = NP*rc_P;
vc_N = NP*vc_P;

end

