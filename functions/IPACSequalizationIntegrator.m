function [time, Xvec, RNvec, BRvec, H_Nvec, Tvec, commandedRates_vec, servoTracking, torques] = IPACSequalizationIntegrator(X0, N, t0, tmax, Gs_B_t0, Gt_B_t0, Gg_B, inertia, gains)
    % Constant inertias
%     Is_B = diag([86, 85, 113]); % kgm^2
%     J_G = diag([0.13, 0.04, 0.03]); % kgm^2 IG + Iw
%     Iws = 0.1; % kgm^2
    Is_B = inertia.Is_B; % kgm^2
    J_G = inertia.J_G; % kgm^2 IG + Iw
    Iws = inertia.Iws; % kgm^2
    
    Js = J_G(1,1);
    Jt = J_G(2,2);
    Jg = J_G(3,3);
%     IGs = Js - Iws;
    
    % Set initial conditions and preallocate vectors
    delta_t = 0.01; %[s]
    
    time = t0:delta_t:tmax;
    Xvec = zeros(3*N+6,length(time));
    RNvec = zeros(6,length(time));
    BRvec = zeros(6,length(time));
    commandedRates_vec = zeros(N*2,length(time));
    H_Nvec = zeros(3,length(time));
    Tvec = zeros(1,length(time));
    servoTracking = zeros(4*N,length(time));
    torques = zeros(2*N+3, length(time));
    
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

        [Lr_B, sigRN, omegaRN_B, sigBR, omegaBR_B] = requiredTorque(time(n), X, N, I_B, Iws, Gs_B, Gt_B, Gg_B, gains);

        [d_OMEGA_desired, d_gamma_desired] = IPACSequalizationCommandedRates(time(n), X, Lr_B, N, Iws, J_G, Gs_B, Gt_B, Gg_B, gains);

        % Calculate control torques
%         K_gamma = 10;
        K_gamma = gains.gamma;
        
        d2_gamma_desired = (d_gamma_desired - commandedRates_vec(1:N,n))/delta_t;
        delta_d_gamma = d_gamma - d_gamma_desired;
        d2_gamma = d2_gamma_desired - K_gamma*delta_d_gamma;
        
        servoTracking(1:2*N,n) = [d_gamma; d_gamma_desired]; % save for servo tracking plots
        
        ug = zeros(N,1);
        us = zeros(N,1);
        for i = 1:N
            omegas = Gs_B(:,i)'*omegaBN_B;
            omegat = Gt_B(:,i)'*omegaBN_B;
            
            ug(i) = Jg*(d2_gamma(i)) - (Js - Jt)*omegas*omegat - Iws*OMEGA(i)*omegat;
            us(i) = Iws*(d_OMEGA_desired(i) + d_gamma(i)*omegat);
        end
        u = [ug; us];
        
        %RK4
        k1 = delta_t.*Xdot(X, u, L_B, N, Gs_B_t0, Gt_B_t0, Gg_B, gamma_t0, inertia);
        k2 = delta_t.*Xdot(X + k1/2, u, L_B, N, Gs_B_t0, Gt_B_t0, Gg_B, gamma_t0, inertia);
        k3 = delta_t.*Xdot(X + k2/2, u, L_B, N, Gs_B_t0, Gt_B_t0, Gg_B, gamma_t0, inertia);
        k4 = delta_t.*Xdot(X + k3, u, L_B, N, Gs_B_t0, Gt_B_t0, Gg_B, gamma_t0, inertia);
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
        torques(:,n) = [us; ug; Lr_B];

    end
    
    sigBN = X(1:3);
%     omegaBN_B = X(4:6);
    [sigRN, omegaRN_R] = missionTracking(time(n+1));
    BN = MRP2C(sigBN);
    RN = MRP2C(sigRN);
    BR = BN*RN';
    omegaRN_B = BR*omegaRN_R;
    RNvec(:,n+1) = [sigRN; omegaRN_B];

end