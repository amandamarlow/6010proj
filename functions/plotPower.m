function [] = plotPower(t, P_desired, P_actual)
%PLOTPOWER Summary of this function goes here
%   Detailed explanation goes here

figure
hold on
plot(t, P_desired, 'k');
plot(t, P_actual,'--r')
title("Power Profile")
ylabel("Power [W]", 'Interpreter', 'Latex')
xlabel("Time [s]")
legend("desired", "actual")
grid on
hold off

end

