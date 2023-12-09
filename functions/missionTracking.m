function [sigRN, omegaRN_R] = missionTracking(t)
    f = 0.03; % [rad/s]
    sigRN = 1/4*[0.1*sin(f*t); 0.2*cos(f*t); -0.3*sin(2*f*t)];
    d_sigRN = 1/4*[0.1*cos(f*t)*f; -0.2*sin(f*t)*f; -0.3*cos(2*f*t)*2*f];
    omegaRN_R = 4*BinvMRP(sigRN)*d_sigRN; 
end

