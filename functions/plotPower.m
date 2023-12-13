function [] = plotPower(t, P_desired, P_actual)
%PLOTPOWER Summary of this function goes here
%   Detailed explanation goes here

figure
hold on
plot(t, P_desired, 'k');
scatter(t, P_actual,'k*')
title("Power Profile")
ylabel("Power [W]", 'Interpreter', 'Latex')
xlabel("Time [s]")
grid on
hold off

end

