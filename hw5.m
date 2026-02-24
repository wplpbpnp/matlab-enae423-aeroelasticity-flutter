% Benjamin Stutzke
% ENAE 423
% Homework 5

R = 0.0461;

K = [(6/5)+R    -1/5    -R;
      -1/5       1/5     0;
       -R         0      R;];

M = [1 0   0;
     0 1/4 0;
     0 0   1/12;];

F = [0 1 0]';
q_static = K\F;

omega_bar = 0:0.00001:2;
for i=1:length(omega_bar)
    A(:, i) = (K-((omega_bar(i).^2).*M))\F;
    A_bar(:, i) = A(:,i)/q_static(1);
    Amp(:, i) = abs(A_bar(:, i));
end

figure(1);
plot(omega_bar, A_bar(1, :));
title('$\bar{A_1}$ for 3DOF System','Interpreter','Latex');
xlabel('$\bar{\Omega}$','Interpreter','Latex');
ylabel('$\bar{A}$','Interpreter','Latex');
ylim([-100 100]);
grid on;

figure(2);
plot(omega_bar, A_bar(2, :));
title('$\bar{A_2}$ for 3DOF System','Interpreter','Latex');
xlabel('$\bar{\Omega}$','Interpreter','Latex');
ylabel('$\bar{A}$','Interpreter','Latex');
ylim([-100 100]);
grid on;

figure(3);
plot(omega_bar, A_bar(3, :));
title('$\bar{A_3}$ for 3DOF System','Interpreter','Latex');
xlabel('$\bar{\Omega}$','Interpreter','Latex');
ylabel('$\bar{A}$','Interpreter','Latex');
ylim([-100 100]);
grid on;

figure(4);
plot(omega_bar, Amp(1, :));
title('$\bar{q_1}$ for 3DOF System','Interpreter','Latex');
xlabel('$\bar{\Omega}$','Interpreter','Latex');
ylabel('$\bar{q}$','Interpreter','Latex');
ylim([0 100]);
grid on;

figure(5);
plot(omega_bar, Amp(2, :));
title('$\bar{q_2}$ for 3DOF System','Interpreter','Latex');
xlabel('$\bar{\Omega}$','Interpreter','Latex');
ylabel('$\bar{q}$','Interpreter','Latex');
ylim([0 100]);
grid on;

figure(6);
plot(omega_bar, Amp(3, :));
title('$\bar{q_3}$ for 3DOF System','Interpreter','Latex');
xlabel('$\bar{\Omega}$','Interpreter','Latex');
ylabel('$\bar{q}$','Interpreter','Latex');
ylim([0 100]);
grid on;
