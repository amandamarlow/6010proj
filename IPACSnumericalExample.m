%% Simulation of VSCMG Null Motion 
% Amanda Marlow 10/20/23

clear
clc
close all
addpath("C:\Users\marlo\MATLAB Drive\6010\RigidBodyKinematics-Matlab\Matlab")

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
[time, Xvec, RNvec, BRvec, H_Nvec, Tvec, commandedRates_vec, torques, servoTracking] = integrator(X0, N, 0, tf, gamma_tf, Gs_B_t0, Gt_B_t0, Gg_B_t0, Is_B);

% Plot
plotAllVSCMG(time, N, Xvec, RNvec, BRvec, H_Nvec, Tvec, commandedRates_vec, torques, servoTracking)

%% RK4 Integrator
function [time, Xvec, RNvec, BRvec, H_Nvec, Tvec, commandedRates_vec, torques, servoTracking] = integrator(X0, N, t0, tmax, gamma_tf, Gs_B_t0, Gt_B_t0, Gg_B, Is_B)
    
    % set all gains
    K = 5;
    P = 15;
    K_gamma = 10;
    ke = 10;
    
    sigRN = X0(1:3);
    omegaRN_B = zeros(3,1);

    % Constant inertias
    J_G = diag([0.13, 0.04, 0.03]); % kgm^2 IG + Iw
    Iws = 0.1; % kgm^2
    
    Js = J_G(1,1);
    Jt = J_G(2,2);
    Jg = J_G(3,3);
    
    % Set initial conditions and preallocate vectors
    delta_t = 0.01; %[s]
    
    time = t0:delta_t:tmax;
    Xvec = zeros(3*N+6,length(time));
    RNvec = zeros(6,length(time));
    BRvec = zeros(6,length(time));
    commandedRates_vec = zeros(N*2,length(time));
    H_Nvec = zeros(3,length(time));
    Tvec = zeros(1,length(time));
    torques = zeros(2*N+3,length(time));
    servoTracking = zeros(4*N,length(time));
    
    gamma_t0 = X0(7:6+N);
    
    % calculate angular momentumand rotational kinetic energy of initial state
    [H_N_t0, T_t0] = getMomentumEnergy(X0, N, Is_B, J_G, Iws, Gs_B_t0, Gt_B_t0, Gg_B, gamma_t0);
    
    % Include initial conditions in output vectors
    Xvec(:,1) = X0;
    H_Nvec(:,1) = H_N_t0;
    Tvec(:,1) = T_t0;
    
    % Define External Torque
    L_B = 0;
    
    X = X0;
    for n = 1:length(time)-1
        
        % extract current state
        sigBN = X(1:3);
        omegaBN_B = X(4:6);
        gamma = X(7:6+N);
        d_gamma = X(7+N:6+2*N);
        OMEGA = X(7+2*N:6+3*N);
        
        Gs_B = zeros(3,N);
        Gt_B = zeros(3,N);
        J_B_sum = 0;
        for i = 1:N
            gs_B = cos(gamma(i) - gamma_t0(i))*Gs_B_t0(:,i) + sin(gamma(i) - gamma_t0(i))*Gt_B_t0(:,i);
            gt_B = -sin(gamma(i) - gamma_t0(i))*Gs_B_t0(:,i) + cos(gamma(i) - gamma_t0(i))*Gt_B_t0(:,i);
            gg_B = Gg_B(:,i);

            Gs_B(:,i) = gs_B;
            Gt_B(:,i) = gt_B;
        
            % assumes all vscmgs have same J_G
            J_B_sum = J_B_sum + Js*(gs_B*gs_B') + Jt*(gt_B*gt_B') + Jg*(gg_B*gg_B');
           
        end
        I_B = Is_B + J_B_sum;
        
        % calculate required torque and attitude error
        [Lr_B, sigRN, omegaRN_B, sigBR, omegaBR_B] = requiredTorque(sigRN, X, N, I_B, Iws, Gs_B, Gt_B, Gg_B, K, P);
        % outer loop servo rate commands
        [d_eta] = commandedRates(X, Lr_B, N, Iws, J_G, Gs_B, Gt_B, Gg_B, gamma_tf, ke);
        
        d_gamma_desired = d_eta(N+1:end);
        d_OMEGA_desired = d_eta(1:N);
        d2_gamma_desired = (d_gamma_desired - commandedRates_vec(1:N,n))/delta_t;
        
        % Calculate control torques
        delta_d_gamma = d_gamma - d_gamma_desired;
        d2_gamma = d2_gamma_desired - K_gamma*delta_d_gamma;
        
        servoTracking(1:2*N,n) = [d_gamma; d_gamma_desired]; % save for servo tracking plots
        
        ug = zeros(N,1);
        us = zeros(N,1);
        for i = 1:N

            omegas = Gs_B(:,i)'*omegaBN_B;
            omegat = Gt_B(:,i)'*omegaBN_B;
%             omegag = gg_B'*omegaBN_B;
            
            % PROBLEM
            ug(i) = Jg*(d2_gamma(i)) - (Js - Jt)*omegas*omegat - Iws*OMEGA(i)*omegat;
            us(i) = Iws*(d_OMEGA_desired(i) + d_gamma(i)*omegat);
        end
        u = [ug; us];
        
        %RK4
        k1 = delta_t.*Xdot(X, u, L_B, N, Gs_B_t0, Gt_B_t0, Gg_B, gamma_t0);
        k2 = delta_t.*Xdot(X + k1/2, u, L_B, N, Gs_B_t0, Gt_B_t0, Gg_B, gamma_t0);
        k3 = delta_t.*Xdot(X + k2/2, u, L_B, N, Gs_B_t0, Gt_B_t0, Gg_B, gamma_t0);
        k4 = delta_t.*Xdot(X + k3, u, L_B, N, Gs_B_t0, Gt_B_t0, Gg_B, gamma_t0);
        X = X + 1/6*(k1 + 2*k2 + 2*k3 + k4);
        
        if norm(X(1:3))>1
            X(1:3) = -X(1:3)/(X(1:3)'*X(1:3));
        end
        
        % calculate angular momentumand rotational kinetic energy
        [H_N, T] = getMomentumEnergy(X, N, Is_B, J_G, Iws, Gs_B_t0, Gt_B_t0, Gg_B, gamma_t0);
        
        % get delta omega dot
        d_X = 1/6*(k1 + 2*k2 + 2*k3 + k4)/delta_t;
        d_OMEGA = d_X(end-(N-1):end);
        servoTracking(2*N+1:end,n) = [d_OMEGA; d_OMEGA_desired];
        
        % store all info in output vectors
        Xvec(:,n+1) = X;
        H_Nvec(:,n+1) = H_N;
        Tvec(:,n+1) = T;
        RNvec(:,n) = [sigRN; omegaRN_B];
        BRvec(:,n) = [sigBR; omegaBR_B];
        commandedRates_vec(:,n+1) = [d_gamma_desired; d_OMEGA_desired];
        torques(:,n+1) = [us; ug; Lr_B];

    end
    
    RNvec(:,n+1) = [sigRN; omegaRN_B];
    BRvec(:,n) = [sigBR; omegaBR_B];

end
%% Rates Function
function [Xdot] = Xdot(X, u, L_B, N, Gs_B_t0, Gt_B_t0, Gg_B_t0, gamma_t0)
    
    Is_B = diag([86, 85, 113]); % kgm^2
    J_G = diag([0.13, 0.04, 0.03]); % kgm^2 IG + Iw
    Iws = 0.1; % kgm^2
    
    sigBN = X(1:3);
    omegaBN_B = X(4:6);
    gamma = X(7:6+N);
    d_gamma = X(7+N:6+2*N);
    OMEGA = X(7+2*N:6+3*N);
    
    ug = u(1:N);
    us = u(N+1:end);
    
    Js = J_G(1,1);
    Jt = J_G(2,2);
    Jg = J_G(3,3);
    
    Gs_B = zeros(3,N);
    Gt_B = zeros(3,N);
    Gg_B = Gg_B_t0;
    
    f_omega_sum = 0;
    f_d_gamma = zeros(N,1);
    f_OMEGA = zeros(N,1);
    J_B_sum = 0;
    for i = 1:N
        gs_B = cos(gamma(i) - gamma_t0(i))*Gs_B_t0(:,i) + sin(gamma(i) - gamma_t0(i))*Gt_B_t0(:,i);
        gt_B = -sin(gamma(i) - gamma_t0(i))*Gs_B_t0(:,i) + cos(gamma(i) - gamma_t0(i))*Gt_B_t0(:,i);
        gg_B = Gg_B_t0(:,i);
        
        Gs_B(:,i) = gs_B;
        Gt_B(:,i) = gt_B;
        
        % assumes all vscmgs have same J_G
        J_B_sum = J_B_sum + Js*(gs_B*gs_B') + Jt*(gt_B*gt_B') + Jg*(gg_B*gg_B');
        
        omegas = gs_B'*omegaBN_B;
        omegat = gt_B'*omegaBN_B;
        omegag = gg_B'*omegaBN_B;
        
        f_omega_sum =  f_omega_sum ...
            + gs_B*(Js*d_gamma(i)*omegat - (Jt-Jg)*omegat*d_gamma(i))...
            + gt_B*((Js*omegas + Iws*OMEGA(i))*d_gamma(i) - (Jt+Jg)*omegas*d_gamma(i) + Iws*OMEGA(i)*omegag)...
            - gg_B*Iws*OMEGA(i)*omegat;
        
        f_d_gamma(i) = ug(i) + (Js - Jt)*omegas*omegat + Iws*OMEGA(i)*omegat;
        f_OMEGA(i) = us(i) - Iws*d_gamma(i)*omegat;
    end
    
    I_B = Is_B + J_B_sum;
    
    B = ((1-sigBN'*sigBN)*eye(3)+2*tilde(sigBN)+2*(sigBN*sigBN'));
    f_sig = 0.25*B*omegaBN_B;    
    f_omega = -tilde(omegaBN_B)*I_B*omegaBN_B + L_B - f_omega_sum;
    f_gamma = d_gamma;
    
    fmat = [f_sig; f_omega; f_gamma; f_d_gamma; f_OMEGA];
    
    M = [eye(3,3), zeros(3,3), zeros(3,N), zeros(3,N), zeros(3,N); ...
        zeros(3,3), I_B, zeros(3,N), Jg*Gg_B, Iws*Gs_B; ...
        zeros(N,3), zeros(N,3), eye(N,N), zeros(N,N), zeros(N,N); ...
        zeros(N,3), Jg*Gg_B', zeros(N,N), Jg*eye(N,N), zeros(N,N); ...
        zeros(N,3), Iws*Gs_B', zeros(N,N), zeros(N,N), Iws*eye(N,N)];
    
    Xdot = M \ fmat;
end

%% get angular momentum and kinetic energy
function [H_N, T] = getMomentumEnergy(X, N, IRW_B, J_G, Iws, Gs_B_t0, Gt_B_t0, Gg_B_t0, gamma_t0)

    % calculate inertial angular momentum and kinetic energy from state
    sigBN = X(1:3);
    omegaBN_B = X(4:6);
    gamma = X(7:6+N);
    d_gamma = X(7+N:6+2*N);
    OMEGA = X(7+2*N:6+3*N);
    
    Js = J_G(1,1);
    Jt = J_G(2,2);
    Jg = J_G(3,3);
    IGs = Js - Iws;

    HB_B = IRW_B*omegaBN_B;

    HGW_B = 0;
    T_sum = 0;
    for i = 1:N
        gs_B = cos(gamma(i) - gamma_t0(i))*Gs_B_t0(:,i) + sin(gamma(i) - gamma_t0(i))*Gt_B_t0(:,i);
        gt_B = -sin(gamma(i) - gamma_t0(i))*Gs_B_t0(:,i) + cos(gamma(i) - gamma_t0(i))*Gt_B_t0(:,i);
        gg_B = Gg_B_t0(:,i);

        omegas = gs_B'*omegaBN_B;
        omegat = gt_B'*omegaBN_B;
        omegag = gg_B'*omegaBN_B;

        HGW_B = HGW_B + Iws*OMEGA(i)*gs_B;

        T_sum = T_sum + Iws*(OMEGA(i) + omegas)^2 + IGs*omegas^2 + Jt*omegat^2 + Jg*(omegag + d_gamma(i))^2;
    end

    H_B = HB_B + HGW_B;
    BN = MRP2C(sigBN);
    NB = BN';
    H_N = NB*H_B;

    T = 0.5*omegaBN_B'*IRW_B*omegaBN_B + 1/2*T_sum;
        
end

%% Outer loop desired gimbal and spin rate function
function [d_eta] = commandedRates(X, Lr_B, N, Iws, J_G, Gs_B, Gt_B, Gg_B, gamma_tf, ke)
    sigBN = X(1:3);
    omegaBN_B = X(4:6);
    gamma = X(7:6+N);
    d_gamma = X(7+N:6+2*N);
    OMEGA = X(7+2*N:6+3*N);

    % get reference state omega RN_B (regulation)
    omegaRN_B = zeros(3,1);
    
    % other stuff
    Js = J_G(1,1);
    Jt = J_G(2,2);
    Jg = J_G(3,3);
    
    D0 = zeros(3,N);
    D1 = zeros(3,N);
    D2 = zeros(3,N);
    D3 = zeros(3,N);
    D4 = zeros(3,N);
    
    for i = 1:N
        omegas = Gs_B(:,i)'*omegaBN_B;
        omegat = Gt_B(:,i)'*omegaBN_B;
        omegag = Gg_B(:,1)'*omegaBN_B;

        D0(:,i) = Gs_B(:,i)*Iws;
        D1(:,i) = Gt_B(:,i)*(Iws*OMEGA(i) + Js/2*omegas) + Js/2*omegat*Gs_B(:,i);
        D2(:,i) = 1/2*Jt*(omegat*Gs_B(:,i) + omegas*Gt_B(:,i));
        D3(:,i) = Jg*(omegat*Gs_B(:,i) - omegas*Gt_B(:,i));
        D4(:,i) = 1/2*(Js - Jt)*(Gs_B(:,i)*Gt_B(:,i)'*omegaRN_B + Gt_B(:,i)*Gs_B(:,i)'*omegaRN_B); %% SHOULD THIS BE IN R COORDINATES
    end
    D = D1 - D2 + D3 + D4;

    % Get desired rates
    mu = 10^(-9);
    h = Js*14.4;
    delta = det(1/(h^2)*(D1*D1'));
    Ws = zeros(1,N);
    for i = 1:N
        Ws(i) = 200*exp(-mu*delta);
    end
    W = diag([Ws, ones(1,N)]);
    Q = [D0, D];
    d_eta = W*Q'*((Q*W*Q')\(-Lr_B));
    d_OMEGA_desired = d_eta(1:N);
    d_gamma_desired = d_eta(N+1:end);
    
    d_eta_c = [d_OMEGA_desired; d_gamma_desired];
    
    % Null Motion -------------------------------------
    a_rw = 0;
    a_cmg = 1;
    A = [a_rw*eye(N) zeros(N); zeros(N) a_cmg*eye(N)];
    
    delta_OMEGA = ones(N,1);
    delta_gamma = gamma - gamma_tf;
    d_eta_n = ke*((Q'/(Q*Q'))*Q - eye(2*N))*A*[delta_OMEGA; delta_gamma];
    
    % final output
    d_eta = d_eta_c + d_eta_n;
end

%% Compute required control torque and attitude / rate errors
function [Lr_B, sigRN, omegaRN_B, sigBR, omegaBR_B] = requiredTorque(sigRN, X, N, I_B, Iws, Gs_B, Gt_B, Gg_B, K, P)
        
    sigBN = X(1:3);
    omegaBN_B = X(4:6);
    gamma = X(7:6+N);
    d_gamma = X(7+N:6+2*N);
    OMEGA = X(7+2*N:6+3*N);
    
    omegaRN_B = zeros(3,1);
    dN_omegaRN_B = zeros(3,1);
    omegaBR_B = omegaBN_B - omegaRN_B;

    BN = MRP2C(sigBN);
    RN = MRP2C(sigRN);
    BR = BN*RN';

    sigBR = C2MRP(BR);
    
    Lr_B = -K*sigBR - P*omegaBR_B + I_B*(dN_omegaRN_B - tilde(omegaBN_B)*omegaRN_B)...
        + tilde(omegaBN_B)*I_B*omegaBN_B;
    for i = 1:N
        omegas = Gs_B(:,i)'*omegaBN_B;
        omegat = Gt_B(:,i)'*omegaBN_B;
        omegag = Gg_B(:,i)'*omegaBN_B;

        Lr_B = Lr_B + Iws*OMEGA(i)*(omegag*Gt_B(:,i) - omegat*Gg_B(:,i));
    end
end