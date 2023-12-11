%% Simulation of VSCMG Null Motion 
% Amanda Marlow 10/20/23

clear
clc
close all
addpath("C:\Users\marlo\MATLAB Drive\6010\RigidBodyKinematics-Matlab\Matlab")
addpath("functions")

% Simulation Parameters
N = 4; % number of VSCMGs
theta = deg2rad(54.75); % angle eof pyramid sides to base

% Inertias
inertia.Is_B = [15053, 3000, -1000; 3000, 6510, 2000; -1000, 2000, 11122];
inertia.Iws = 0.7;
inertia.J_G = diag([inertia.Iws+0.1, 0.4+0.1, 0.4+0.1]);

% Gains
gains.gamma = 2;
gains.k2 = 2*10^-3;
gains.k3 = 2*10^-3;
% gains.w1 = 1*10^-4; % weight matrix (commandedRates)
gains.w1 = 100; % weight matrix (commandedRates)
gains.w2 = 1; % weight matrix (commandedRates)
gains.K = 5; % on sigBR (requiredTorque)
gains.P = 15; % on omegaBR (requiredTorque)
% gains.K = 50; % on sigBR (requiredTorque) % try not to hardcode use performance measures
% gains.P = 20; % on omegaBR (requiredTorque)

% % Old
% Is_B = diag([86, 85, 113]); % kgm^2
% J_G = diag([0.13, 0.04, 0.03]); % kgm^2 IG + Iw

% gimbal axis matrix
gg1_B_t0 = [cos(theta); 0; sin(theta)];
gg2_B_t0 = [-cos(theta); 0; sin(theta)];
gg3_B_t0 = [0; cos(theta); sin(theta)];
gg4_B_t0 = [0; -cos(theta); sin(theta)];
Gg_B_t0 = [gg1_B_t0, gg2_B_t0, gg3_B_t0, gg4_B_t0];

% spin axis matrix
gs1_B_t0 = [0; 1; 0];
gs2_B_t0 = [0; -1; 0];
gs3_B_t0 = [1; 0; 0];
gs4_B_t0 = [-1; 0; 0];
Gs_B_t0 = [gs1_B_t0, gs2_B_t0, gs3_B_t0, gs4_B_t0];

% transverse axis matrix
gt1_B_t0 = cross(gg1_B_t0, gs1_B_t0);
gt2_B_t0 = cross(gg2_B_t0, gs2_B_t0);
gt3_B_t0 = cross(gg3_B_t0, gs3_B_t0);
gt4_B_t0 = cross(gg4_B_t0, gs4_B_t0);
Gt_B_t0 = [gt1_B_t0, gt2_B_t0, gt3_B_t0, gt4_B_t0];

% Initial State
% sigBN_t0 = [0.1; 0.2; 0.3];
% omegaBN_B_t0 = [0.01; -0.01; 0.005]; % [rad/s]
% gamma_t0 = [0; 0; 90; -90]*pi/180; % rad
% d_gamma_t0 = 0*ones(N,1);
% OMEGA_t0 = 14.4*ones(N,1); % rad/s wheel speed
% X0 = [sigBN_t0; omegaBN_B_t0; gamma_t0; d_gamma_t0; OMEGA_t0];

sigBN_t0 = [0; 0; 0];
omegaBN_B_t0 = [0; 0; 0]; % [rad/s]
gamma_t0 = [pi/2; -pi/2; -pi/2; pi/2]; % rad
d_gamma_t0 = [0; 0; 0; 0];
% OMEGA_t0 = [20000; 17500; 15000; 12500]/2/pi/60; % rad/s wheel speed
OMEGA_t0 = [40; 40; 40; 40]; % rad/s wheel speed
X0 = [sigBN_t0; omegaBN_B_t0; gamma_t0; d_gamma_t0; OMEGA_t0];


% Call Integrator
gamma_tf = [-45, 45, -45, 45]'*pi/180; % rad
% tf = 1000; % s
tf = 100;
% [time, Xvec, RNvec, BRvec, H_Nvec, Tvec, commandedRates_vec, servoTracking, torques] = integrator(X0, N, 0, tf, Gs_B_t0, Gt_B_t0, Gg_B_t0, inertia, gains);
[time, Xvec, RNvec, BRvec, H_Nvec, Tvec, commandedRates_vec, servoTracking, torques] = IPACSequalizationIntegrator(X0, N, 0, tf, Gs_B_t0, Gt_B_t0, Gg_B_t0, inertia, gains);

% % test orbit
% t = linspace(0,2*pi/(14.577788549/24/60^2*2*pi));
% rc_N = zeros(3,length(t));
% eclipse = zeros(1,length(t));
% for i = 1:length(t)
%     [rc_N(:,i), ~] = orbitPropogator(t(i));
%     [~, ~, eclipse(i)] = missionTracking(t(i));
% end
% figure
% plotOrbit(rc_N(:,eclipse==0), "testing orbit")
% plotOrbit(rc_N(:,eclipse==1), "testing orbit")

% Plot
plotAllVSCMG(time, N, Xvec, RNvec, BRvec, H_Nvec, Tvec, commandedRates_vec, torques, servoTracking)