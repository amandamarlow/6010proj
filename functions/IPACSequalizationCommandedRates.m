function [d_OMEGA_desired, d_gamma_desired] = IPACSequalizationCommandedRates(t, X, Lr_B, N, Iws, J_G, Gs_B, Gt_B, Gg_B, gains)
    sigBN = X(1:3);
    omegaBN_B = X(4:6);
    gamma = X(7:6+N);
    d_gamma = X(7+N:6+2*N);
    OMEGA = X(7+2*N:6+3*N);

    % get reference state omega RN_B
    [sigRN, omegaRN_R, eclipse] = missionTracking(t);
    dt_omegaRN = 0.001;
    [sigRN_t2, omegaRN_R_t2, ~] = missionTracking(t + dt_omegaRN);

    BN = MRP2C(sigBN);
    RN = MRP2C(sigRN);
    BR = BN*RN';

    RN2 = MRP2C(sigRN_t2);
    BR2 = BN*RN2';
    omegaRN_B_t2 = BR2*omegaRN_R_t2;

    sigBR = C2MRP(BR);

    omegaRN_B = BR*omegaRN_R;
    omegaBR_B = omegaBN_B - omegaRN_B;
    dN_omegaRN_B = (omegaRN_B_t2 - omegaRN_B) / dt_omegaRN;
    
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
%     mu = 10^(-9);
%     h = Js*14.4;
%     h = Js*40;
%     delta = det(1/(h^2)*(D1*D1'));
%     Ws = zeros(1,N);
%     for i = 1:N
%         Ws(i) = 200*exp(-mu*delta);
%     end
%     W = diag([Ws, ones(1,N)]);
    Ws = gains.w1*exp(-gains.w2*cond(D));
    W = diag([Ws*ones(1,N), ones(1,N)]);
%     W = diag([ones(1,N), Ws*ones(1,N)]);

    % Look at singular value decomposition

    P_eclipse = -680; %[W] wuthdraw during eclipse
    P_sun = 1000; %[W]
    if eclipse == 1
        P = P_eclipse;
    else
        P = P_sun;
    end
    
    Q = [D, D0; zeros(1,N), OMEGA'*Iws];
    Lrp = [-Lr_B; P];
    
    OMEGAavg = sum(OMEGA)/N;
%     OMEGAavg = 40;
    OMEGAe = OMEGA - OMEGAavg*ones(N,1);
    R = [zeros(1,N), gains.k3*OMEGAe'];
    u = W * ((Q'*(Q*W*Q')^-1)*(Lrp + Q*W*R') - R');
    
    d_OMEGA_desired = u(N+1:end);
    d_gamma_desired = u(1:N);
end