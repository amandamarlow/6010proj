function [] = plotAllVSCMG(time, N, Xvec, RNvec, BRvec, H_Nvec, Tvec, commandedRates_vec, torques, servoTracking)
%PLOTALLVSCMG Summary of this function goes here
%   Detailed explanation goes here

% Body Attitude and Rates
figure
subplot(2,1,1)
hold on
plot(time, Xvec(1:3,:))
plot(time, RNvec(1:3, :), '--');
sgtitle("Body Attitude and Rates")
ylabel("$\sigma_{BN}$", 'Interpreter', 'Latex')
hold off

subplot(2,1,2)
hold on
plot(time, Xvec(4:6,:))
plot(time, RNvec(4:6, :), '--')
ylabel("$^B \omega_{BN}$ [rad/s]", 'Interpreter', 'Latex')
hold off

% Integrals of Motion
figure
subplot(2,1,1)
plot(time, H_Nvec)
sgtitle("Integrals of Motion")
ylabel("$^N H$ [Nms]", 'Interpreter', 'Latex')

subplot(2,1,2)
plot(time, Tvec)
ylim([min(Tvec)-0.1, max(Tvec)+0.1])
ylabel("$T$ [J]", 'Interpreter', 'Latex')
xlabel("time [s]")

% VSCMG States
figure
subplot(2,1,1)
plot(time, Xvec(7:6+N, :))
ylabel("$\gamma$ [rad]", 'Interpreter', 'Latex')
sgtitle("VSCMG States")

subplot(2,1,2)
plot(time, Xvec(7+2*N:6+3*N, :))
xlabel("time [s]")
ylabel("$\Omega$ [rad/s]", 'Interpreter', 'Latex')

% Attitude Error
figure
subplot(2,1,1)
semilogy(time, vecnorm(BRvec(1:3,:)))
sgtitle("Attitude Error")
ylabel("$|\sigma_{B/R}|$", 'Interpreter', 'Latex')

subplot(2,1,2)
semilogy(time, vecnorm(BRvec(4:6,:)))
ylabel("$|\omega_{B/R}|$ [rad/s]", 'Interpreter', 'Latex') %% why does hint say omega_B/BR
xlabel("Time [s]")

% Control Torques
figure
subplot(3,1,1)
plot(time, torques(1:N,:))
sgtitle("Control Torques")
ylabel("$u_s$ [Nm]", 'Interpreter', 'Latex')

subplot(3,1,2)
plot(time, torques(N+1:2*N,:))
ylabel("$u_g$ [Nm]", 'Interpreter', 'Latex')

subplot(3,1,3)
plot(time, torques(2*N+1:end,:))
xlabel("time [s]")
ylabel("$^B L_r$ [Nm]", 'Interpreter', 'Latex')

% Servo Tracking
figure
subplot(2,2,1)
hold on
plot(time, servoTracking(1:N, :))
plot(time, servoTracking(N+1:2*N, :), '--')
sgtitle("Servo Tracking")
ylabel("$\dot{\gamma}$ [rad/s]", 'Interpreter', 'Latex')
hold off

subplot(2,2,2)
hold on
plot(time, servoTracking(2*N+1:3*N, :))
plot(time, servoTracking(3*N+1:end, :), '--')
ylabel("$\dot{\Omega}$ [rad/s]", 'Interpreter', 'Latex')
hold off

subplot(2,2,3)
semilogy(time, abs(servoTracking(1:N, :) - servoTracking(N+1:2*N, :)))
ylabel("$|\Delta \dot{\gamma}|$ [rad/s]", 'Interpreter', 'Latex')
xlabel("time [s]")

subplot(2,2,4)
semilogy(time, abs(servoTracking(2*N+1:3*N, :) - servoTracking(3*N+1:end, :)))
ylabel("$|\Delta \dot{\Omega}|$ [rad/s]", 'Interpreter', 'Latex')
xlabel("time [s]")
end

