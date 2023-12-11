function [] = plotPower(t, P_desired, P_actual)
%PLOTPOWER Summary of this function goes here
%   Detailed explanation goes here

figure
hold on
plot(t, P_actual)
plot(t, P_desired, '--');
title("Power Profile")
ylabel("Power [W]", 'Interpreter', 'Latex')
xlabel("Time [s]")
grid on
hold off

end

