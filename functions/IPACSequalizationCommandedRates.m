function [d_OMEGA_desired, d_gamma_desired, P, condition] = IPACSequalizationCommandedRates(t, X, Lr_B, N, Iws, J_G, Gs_B, Gt_B, Gg_B, gains, stored_power)
    sigBN = X(1:3);
    omegaBN_B = X(4:6);
    gamma = X(7:6+N);
    d_gamma = X(7+N:6+2*N);
    OMEGA = X(7+2*N:6+3*N);

    % get reference state omega RN_B
    [sigRN, omegaRN_R, P] = missionTracking(t);
    if stored_power >= 1500*60^2; % [Ws]
       P = 0; 
    end
%     dt_omegaRN = 0.001;
%     [sigRN_t2, omegaRN_R_t2, ~] = missionTracking(t + dt_omegaRN);

    BN = MRP2C(sigBN);
    RN = MRP2C(sigRN);
    BR = BN*RN';

%     RN2 = MRP2C(sigRN_t2);
%     BR2 = BN*RN2';
%     omegaRN_B_t2 = BR2*omegaRN_R_t2;

%     sigBR = C2MRP(BR);

    omegaRN_B = BR*omegaRN_R;
%     omegaBR_B = omegaBN_B - omegaRN_B;
%     dN_omegaRN_B = (omegaRN_B_t2 - omegaRN_B) / dt_omegaRN;
    
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
    condition = cond(D);
    % Get desired rates
%     mu = 10^(-9);
%     h = Js*14.4;
%     h = Js*40;
%     delta = det(1/(h^2)*(D1*D1'));
%     Ws = zeros(1,N);
%     for i = 1:N
%         Ws(i) = 200*exp(-mu*delta);
%     end
%     W = diag([Ws, ones(1,N)]);
    Ws = gains.w1*exp(-gains.w2*condition);
    W = diag([Ws*ones(1,N), ones(1,N)]);
    
    
%     % wheel speed equalization with "gradient method" ---------------------
%     Lrp = [-Lr_B; P];
%     Q = [D, D0; zeros(1,N), OMEGA'*Iws];
%     OMEGAavg = sum(OMEGA)/N;
%     OMEGAe = OMEGA - OMEGAavg*ones(N,1);
%     R = [zeros(1,N), gains.k3*OMEGAe'];
%     u = W * ((Q'*(Q*W*Q')^-1)*(Lrp + Q*W*R') - R');
%     % ---------------------------------------------------------------------
    
% wheel speed equalization with extra constraint --------------------------
    OMEGAavg = sum(OMEGA)/N;
    OMEGAe = OMEGA - OMEGAavg*ones(N,1);
    Jw1 = OMEGAe'*OMEGAe/2;
    delta_Jw1 = OMEGAe';

    Lrp = [-Lr_B; P; -gains.k2*Jw1];
    Q = [D, D0; zeros(1,N), OMEGA'*Iws; zeros(1,N), delta_Jw1];
    
    if cond(Q) > 70 %rank(Q) < 5
        u = W^(1/2)*pinv(Q*W^(1/2), 1e-12)*Lrp;
    else
        u = W*Q'*((Q*W*Q')^-1)*(Lrp);
    end
    
% -------------------------------------------------------------------------
    
    d_OMEGA_desired = u(N+1:end);
    d_gamma_desired = u(1:N);
end