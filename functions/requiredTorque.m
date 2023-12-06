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

