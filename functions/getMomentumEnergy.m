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
