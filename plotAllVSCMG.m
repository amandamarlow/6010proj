function [] = plotAllVSCMG(time, N, Xvec, RNvec, BRvec, H_Nvec, Tvec, commandedRates_vec, torques, servoTracking)
%PLOTALLVSCMG Summary of this function goes here
%   Detailed explanation goes here

% Body Attitude and Rates
figure
subplot(2,1,1)
hold on
plot(time, Xvec(1:3,:))
set(gca,'ColorOrderIndex',1)
plot(time, RNvec(1:3, :), '--');
sgtitle("Body Attitude and Rates")
ylabel("$\sigma_{BN}$", 'Interpreter', 'Latex')
legend('$\sigma_1$', '$\sigma_2$', '$\sigma_3$', 'Location', 'best', 'Interpreter', 'Latex')
grid on
hold off

subplot(2,1,2)
hold on
plot(time, Xvec(4:6,:))
set(gca,'ColorOrderIndex',1)
plot(time, RNvec(4:6, :), '--')
ylabel("$^B \omega_{BN}$ [rad/s]", 'Interpreter', 'Latex')
legend('$\omega_1$', '$\omega_2$', '$\omega_3$', 'Location', 'best', 'Interpreter', 'Latex')
grid on
hold off

% Integrals of Motion
figure
subplot(2,1,1)
plot(time, H_Nvec)
sgtitle("Integrals of Motion")
ylabel("$^N H$ [Nms]", 'Interpreter', 'Latex')
legend('$H_1$', '$H_2$', '$H_3$', 'Location', 'best', 'Interpreter', 'Latex')
grid on

subplot(2,1,2)
plot(time, Tvec)
ylim([min(Tvec)-0.1, max(Tvec)+0.1])
ylabel("$T$ [J]", 'Interpreter', 'Latex')
xlabel("time [s]")
grid on

% VSCMG States
figure
subplot(2,1,1)
plot(time, Xvec(7:6+N, :))
ylabel("$\gamma$ [rad]", 'Interpreter', 'Latex')
sgtitle("VSCMG States")
legend('$\gamma_1$', '$\gamma_2$', '$\gamma_3$', '$\gamma_4$', 'Location', 'best', 'Interpreter', 'Latex')
grid on

subplot(2,1,2)
plot(time, Xvec(7+2*N:6+3*N, :))
xlabel("time [s]")
ylabel("$\Omega$ [rad/s]", 'Interpreter', 'Latex')
legend('$\Omega_1$', '$\Omega_2$', '$\Omega_3$', '$\Omega_4$', 'Location', 'best', 'Interpreter', 'Latex')
grid on

% Attitude Error
figure
subplot(2,1,1)
% semilogy(time, vecnorm(BRvec(1:3,:)))
plot(time, BRvec(1:3,:))
sgtitle("Attitude Error")
legend('$\sigma_1$', '$\sigma_2$', '$\sigma_3$', 'Location', 'best', 'Interpreter', 'Latex')
grid on
ylabel("$|\sigma_{B/R}|$", 'Interpreter', 'Latex')

subplot(2,1,2)
% semilogy(time, vecnorm(BRvec(4:6,:)))
plot(time, BRvec(4:6,:))
ylabel("$|\omega_{B/R}|$ [rad/s]", 'Interpreter', 'Latex') %% why does hint say omega_B/BR
legend('$\omega_1$', '$\omega_2$', '$\omega_3$', 'Location', 'best', 'Interpreter', 'Latex')
xlabel("Time [s]")
grid on

% Control Torques
figure
subplot(3,1,1)
plot(time, torques(1:N,:))
sgtitle("Control Torques")
ylabel("$u_s$ [Nm]", 'Interpreter', 'Latex')
legend('$u_s1$', '$u_s2$', '$u_s3$', '$u_s4$', 'Location', 'best', 'Interpreter', 'Latex')
grid on

subplot(3,1,2)
plot(time, torques(N+1:2*N,:))
ylabel("$u_g$ [Nm]", 'Interpreter', 'Latex')
legend('$u_g1$', '$u_g2$', '$u_g3$', '$u_g4$', 'Location', 'best', 'Interpreter', 'Latex')
grid on

subplot(3,1,3)
plot(time, torques(2*N+1:end,:))
xlabel("time [s]")
ylabel("$^B L_r$ [Nm]", 'Interpreter', 'Latex')
legend('$L_r1$', '$L_r2$', '$L_r3$', 'Location', 'best', 'Interpreter', 'Latex')
grid on

% Servo Tracking
figure
subplot(2,2,1)
hold on
plot(time, servoTracking(1:N, :))
set(gca,'ColorOrderIndex',1)
plot(time, servoTracking(N+1:2*N, :), '--')
legend('$\dot{\gamma}_1$', '$\dot{\gamma}_2$', '$\dot{\gamma}_3$', '$\dot{\gamma}_4$', 'Location', 'best', 'Interpreter', 'Latex')
sgtitle("Servo Tracking")
ylabel("$\dot{\gamma}$ [rad/s]", 'Interpreter', 'Latex')
grid on
hold off

subplot(2,2,2)
hold on
plot(time, servoTracking(2*N+1:3*N, :))
set(gca,'ColorOrderIndex',1)
plot(time, servoTracking(3*N+1:end, :), '--')
ylabel("$\dot{\Omega}$ [rad/s]", 'Interpreter', 'Latex')
legend('$\dot{\Omega}_1$', '$\dot{\Omega}_2$', '$\dot{\Omega}_3$', '$\dot{\Omega}_4$', 'Location', 'best', 'Interpreter', 'Latex')
grid on
hold off

subplot(2,2,3)
semilogy(time, abs(servoTracking(1:N, :) - servoTracking(N+1:2*N, :)))
ylabel("$|\Delta \dot{\gamma}|$ [rad/s]", 'Interpreter', 'Latex')
legend('$\Delta\dot{\gamma}_1$', '$\Delta\dot{\gamma}_2$', '$\Delta\dot{\gamma}_3$', '$\Delta\dot{\gamma}_4$', 'Location', 'best', 'Interpreter', 'Latex')
grid on
xlabel("time [s]")

subplot(2,2,4)
semilogy(time, abs(servoTracking(2*N+1:3*N, :) - servoTracking(3*N+1:end, :)))
legend('$\Delta\dot{\Omega}_1$', '$\Delta\dot{\Omega}_2$', '$\Delta\dot{\Omega}_3$', '$\Delta\dot{\Omega}_4$', 'Location', 'best', 'Interpreter', 'Latex')
ylabel("$|\Delta \dot{\Omega}|$ [rad/s]", 'Interpreter', 'Latex')
xlabel("time [s]")
grid on
end

