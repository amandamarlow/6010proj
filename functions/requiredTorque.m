function [Lr_B, sigRN, omegaRN_B, sigBR, omegaBR_B] = requiredTorque(t, X, N, I_B, Iws, Gs_B, Gt_B, Gg_B, gains)
        
%     K = 5;
%     P = 15;
        
    sigBN = X(1:3);
    omegaBN_B = X(4:6);
    gamma = X(7:6+N);
    d_gamma = X(7+N:6+2*N);
    OMEGA = X(7+2*N:6+3*N);

    [sigRN, omegaRN_R] = missionTracking(t);
    dt_omegaRN = 0.0001;
    [sigRN_t2, omegaRN_R_t2] = missionTracking(t + dt_omegaRN);

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
    
    Lr_B = -gains.K*sigBR - gains.P*omegaBR_B + I_B*(dN_omegaRN_B - tilde(omegaBN_B)*omegaRN_B)...
        + tilde(omegaBN_B)*I_B*omegaBN_B;
    for i = 1:N
        omegas = Gs_B(:,i)'*omegaBN_B;
        omegat = Gt_B(:,i)'*omegaBN_B;
        omegag = Gg_B(:,i)'*omegaBN_B;

        Lr_B = Lr_B + Iws*OMEGA(i)*(omegag*Gt_B(:,i) - omegat*Gg_B(:,i));
    end
end

