function [d_eta] = IPACScommandedRates(t,X, Lr_B, N, Iws, J_G, Gs_B, Gt_B, Gg_B, gamma_tf, ke)

    w1 = 1*10^-4;
    w2 = 1;
%     P = 680; % Required Power
    P = 10*sin(t); 
    
%     Lrp = [-Lr_B; P];
    Lrp = [P; -Lr_B];

    sigBN = X(1:3);
    omegaBN_B = X(4:6);
    gamma = X(7:6+N);
    d_gamma = X(7+N:6+2*N);
    OMEGA = X(7+2*N:6+3*N);

    % get reference state omega RN_B
    [sigRN, omegaRN_R] = RefState(t);
    dt_omegaRN = 0.0001;
    [sigRN_t2, omegaRN_R_t2] = RefState(t + dt_omegaRN);

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
    
%     D0 = zeros(3,N);
    D1 = zeros(3,N);
    D2 = zeros(3,N);
    D3 = zeros(3,N);
    D4 = zeros(3,N);
    
    for i = 1:N
        omegas = Gs_B(:,i)'*omegaBN_B;
        omegat = Gt_B(:,i)'*omegaBN_B;
        omegag = Gg_B(:,1)'*omegaBN_B;

%         D0(:,i) = Gs_B(:,i)*Iws;
        D1(:,i) = Gt_B(:,i)*(Iws*OMEGA(i) + Js/2*omegas) + Js/2*omegat*Gs_B(:,i);
        D2(:,i) = 1/2*Jt*(omegat*Gs_B(:,i) + omegas*Gt_B(:,i));
        D3(:,i) = Jg*(omegat*Gs_B(:,i) - omegas*Gt_B(:,i));
        D4(:,i) = 1/2*(Js - Jt)*(Gs_B(:,i)*Gt_B(:,i)'*omegaRN_B + Gt_B(:,i)*Gs_B(:,i)'*omegaRN_B); %% SHOULD THIS BE IN R COORDINATES
    end
    D = D1 - D2 + D3 + D4;
    D0 = Gs_B*Iws;
    
    if rank([D0, D]) < 3
        disp(rank([D0, D]))
    end
    if rank(D) < 3
        disp(rank(D))
    end
    if rank(D0) < 3
        disp(rank(D0))
    end

    % Get desired rates
    mu = 10^(-9);
    h = Js*14.4;
    delta = det(1/(h^2)*(D1*D1')); %%%%%%%%%%Check that this is the condition number
    S = svd(D);
    condition = S(1)/S(3);
    Ws = zeros(1,N);
    for i = 1:N
%         Ws(i) = w1*exp(-w2*delta);
        Ws(i) = w1*exp(-w2*condition);
    end
    W = diag([Ws, ones(1,N)]);
    Q = [D0, D; OMEGA'*diag(Iws*ones(1,N)), zeros(1,N)];
    d_eta = W*Q'*((Q*W*Q')\(Lrp));

%     disp(rank([D0, D]))

end

