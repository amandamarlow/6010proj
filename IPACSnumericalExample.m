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
omegaBN_B_t0 = [0; 0; 0]; % [rad/s]
% d_omegaBN_B_t0 = [0; 0; 0]; % given in paper, but I'm not using
sigBN_t0 = [0;0;0];
gamma_t0 = [1, -1, -1, 1]'*pi/2; % rad
d_gamma_t0 = 0*ones(N,1);
% WHEEL SPEEDS NOT GIVEN (appear all different in graph
OMEGA_t0 = 14.4*ones(N,1); % rad/s wheel speed
Is_B = [15053, 3000, -1000; 3000, 6510, 2000; -1000, 2000, 11122];
% NEED TO INCORPORATE THEIR Iw and Ig

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
% omegaBN_B_t0 = [0; 0; 0]; % [rad/s]
% gamma_t0 = [0, 0, 90, -90]*pi/180; % rad
% d_gamma_t0 = 0*ones(1,N);
% OMEGA_t0 = 14.4*ones(1,N); % rad/s wheel speed
X0 = [sigBN_t0; omegaBN_B_t0; gamma_t0; d_gamma_t0; OMEGA_t0];

% Call Integrator
gamma_tf = [-45, 45, -45, 45]'*pi/180; % rad
tf = 300; % s
[time, Xvec, RNvec, BRvec, H_Nvec, Tvec, commandedRates_vec, torques, servoTracking] = IPACSintegrator(X0, N, 0, tf, gamma_tf, Gs_B_t0, Gt_B_t0, Gg_B_t0, Is_B);

% Plot
plotAllVSCMG(time, N, Xvec, RNvec, BRvec, H_Nvec, Tvec, commandedRates_vec, torques, servoTracking)