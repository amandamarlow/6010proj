function [Xdot] = Xdot(X, u, L_B, N, Gs_B_t0, Gt_B_t0, Gg_B_t0, gamma_t0, inertia)
    
%     Is_B = diag([86, 85, 113]); % kgm^2
%     J_G = diag([0.13, 0.04, 0.03]); % kgm^2 IG + Iw
%     Iws = 0.1; % kgm^2
    Is_B = inertia.Is_B; % kgm^2
    J_G = inertia.J_G; % kgm^2 IG + Iw
    Iws = inertia.Iws; % kgm^2
    
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
