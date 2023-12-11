function [r_vec, v_vec] = fg(E, E0, r0_vec, v0_vec, t1_t0, mu)
    % r0_V and v0_V are the initial position and velocity vectors. They
    % must be given in a consistent frame and the answer will be in the
    % same frame
    
    r0 = norm(r0_vec);
    v0 = norm(v0_vec);
    
    epsilon = 1/2*v0^2 - mu/r0;
    a = -mu/2/epsilon;
    delta_E = E - E0;
    f = 1 - a/r0*(1-cos(delta_E));
    g = t1_t0 - sqrt(a^3/mu)*(delta_E - sin(delta_E));
    
    r_vec = f*r0_vec + g*v0_vec;
    
    r = norm(r_vec);
    fdot = -sin(delta_E)*sqrt(mu*a)/r0/r;
    gdot = 1-a/r*(1-cos(delta_E));
    
    v_vec = fdot*r0_vec + gdot*v0_vec;
end
