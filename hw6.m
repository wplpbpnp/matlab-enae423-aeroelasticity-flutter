% Benjamin Stutzke
% ENAE 423
% Homework 6

%% Problem 1
Pbar = 0:0.01:6;
lambda1 = zeros(length(Pbar), 1);
lambda2 = lambda1;

for i=1:length(Pbar)
    K_eff = [3-Pbar(i) Pbar(i)-2;
             -2       2  ];
    M = [3 1; 1 1];

    [X, D] = eig(K_eff, M);
    lambda1(i) = D(1, 1);
    lambda2(i) = D(2, 2);

    lambda1(i) = sqrt(lambda1(i)*-1);
    lambda2(i) = sqrt(lambda2(i)*-1);

    %A = M(1,1)*M(2,2)-M(1,2)*M(2,1);
    %B = K_eff(1,1)*M(2,2)+K_eff(2,2)*M(1,1)-K_eff(1,2)*M(2,1)-K_eff(2,1)*M(1,2);
    %C = K_eff(1,1)*K_eff(2,2)-K_eff(1,2)*K_eff(2,1);

    %lambda1(i) = (-B+sqrt(B.^2-4*A*C))/2/A;
    %lambda2(i) = (-B-sqrt(B.^2-4*A*C))/2/A;
end

figure(1);
plot(Pbar, real(lambda1));
title('$Re\{\lambda_1$\} vs $\bar{P}$','Interpreter','Latex');
xlabel('$\bar{P}$','Interpreter','Latex');
ylabel('$Re\{\lambda_1\}$','Interpreter','Latex');
grid on;

figure(2);
plot(Pbar, real(lambda2));
title('$Re\{\lambda_2$\} vs $\bar{P}$','Interpreter','Latex');
xlabel('$\bar{P}$','Interpreter','Latex');
ylabel('$Re\{\lambda_2\}$','Interpreter','Latex');
grid on;

omega1 = sqrt(lambda1);
omega2 = sqrt(lambda2);

figure(3);
hold on
plot(Pbar, abs(imag(omega1)));
plot(Pbar, abs(imag(omega2)));
title('$\omega_1$ and $\omega_2$ vs $\bar{P}$','Interpreter','Latex');
xlabel('$\bar{P}$','Interpreter','Latex');
ylabel('$\omega$','Interpreter','Latex');
legend('$\omega_1$', '$\omega_2$', 'Interpreter', 'Latex');
grid on;

%% Problem 2
m = 8.7;
Sa = 5.4;
Ia = 84;
k1 = 6.3e4;
k2 = 1.3e6;
c = 5;
e = 0.19;
S = 150;
dCLda = 1.87*pi;

rho = 0.002377;

% Case 1
fprintf("Case 1: \n");
% a
q_d = k2/(e*S*dCLda)
U_d = sqrt(2*q_d/rho)

% b
q = 0:0.1:500;

lambda1 = zeros(length(q), 1);
lambda2 = lambda1;
M = [m -Sa; -Sa Ia];

for i=1:length(q)
    
    K_eff = [k1 -q(i)*S*dCLda; 0 k2 - e*q(i)*S*dCLda];

    [X, D] = eig(K_eff, M);
    lambda1(i) = D(1, 1);
    lambda2(i) = D(2, 2);

    lambda1(i) = sqrt(lambda1(i)*-1);
    lambda2(i) = sqrt(lambda2(i)*-1);
end

figure(4);
plot(q, real(lambda1));
title('$Re\{\lambda_1$\} vs $q$ for Case 1','Interpreter','Latex');
xlabel('$q$','Interpreter','Latex');
ylabel('$Re\{\lambda_1\}$','Interpreter','Latex');
grid on;

figure(5);
plot(q, real(lambda2));
title('$Re\{\lambda_2$\} vs $q$ for Case 1','Interpreter','Latex');
xlabel('$q$','Interpreter','Latex');
ylabel('$Re\{\lambda_2\}$','Interpreter','Latex');
grid on;

omega1 = sqrt(lambda1);
omega2 = sqrt(lambda2);

figure(6);
hold on
plot(q, abs(imag(omega1)));
plot(q, abs(imag(omega2)));
title('$\omega_1$ and $\omega_2$ vs $q$ for Case 1','Interpreter','Latex');
xlabel('$q$','Interpreter','Latex');
ylabel('$\omega$','Interpreter','Latex');
legend('$\omega_1$', '$\omega_2$', 'Interpreter', 'Latex');
grid on;

[maximum, index] = max(abs(imag(omega1)));
q_f = q(index)
U_f = sqrt(2*q_f/rho)

% Case 2
fprintf("Case 2: \n");
e = 0.38;
% a
q_d = k2/(e*S*dCLda)
U_d = sqrt(2*q_d/rho)

% b
q = 0:0.1:500;

lambda1 = zeros(length(q), 1);
lambda2 = lambda1;
M = [m -Sa; -Sa Ia];

for i=1:length(q)
    
    K_eff = [k1 -q(i)*S*dCLda; 0 k2 - e*q(i)*S*dCLda];

    [X, D] = eig(K_eff, M);
    lambda1(i) = D(1, 1);
    lambda2(i) = D(2, 2);

    lambda1(i) = sqrt(lambda1(i)*-1);
    lambda2(i) = sqrt(lambda2(i)*-1);
end

figure(7);
plot(q, real(lambda1));
title('$Re\{\lambda_1$\} vs $q$ for Case 2','Interpreter','Latex');
xlabel('$q$','Interpreter','Latex');
ylabel('$Re\{\lambda_1\}$','Interpreter','Latex');
grid on;

figure(8);
plot(q, real(lambda2));
title('$Re\{\lambda_2$\} vs $q$ for Case 2','Interpreter','Latex');
xlabel('$q$','Interpreter','Latex');
ylabel('$Re\{\lambda_2\}$','Interpreter','Latex');
grid on;

omega1 = sqrt(lambda1);
omega2 = sqrt(lambda2);

figure(9);
hold on
plot(q, abs(imag(omega1)));
plot(q, abs(imag(omega2)));
title('$\omega_1$ and $\omega_2$ vs $q$ for Case 2','Interpreter','Latex');
xlabel('$q$','Interpreter','Latex');
ylabel('$\omega$','Interpreter','Latex');
legend('$\omega_1$', '$\omega_2$', 'Interpreter', 'Latex');
grid on;

[maximum, index] = max(abs(imag(omega1)));
q_f = q(index)
U_f = sqrt(2*q_f/rho)

% Case 3
fprintf("Case 3: \n");
e = 0.19;
Sa = 1.5;

% a
q_d = k2/(e*S*dCLda)
U_d = sqrt(2*q_d/rho)

% b
q = 0:0.1:1000;

lambda1 = zeros(length(q), 1);
lambda2 = lambda1;
M = [m -Sa; -Sa Ia];

for i=1:length(q)
    
    K_eff = [k1 -q(i)*S*dCLda; 0 k2 - e*q(i)*S*dCLda];

    [X, D] = eig(K_eff, M);
    lambda1(i) = D(1, 1);
    lambda2(i) = D(2, 2);

    lambda1(i) = sqrt(lambda1(i)*-1);
    lambda2(i) = sqrt(lambda2(i)*-1);
end

figure(10);
plot(q, real(lambda1));
title('$Re\{\lambda_1$\} vs $q$ for Case 3','Interpreter','Latex');
xlabel('$q$','Interpreter','Latex');
ylabel('$Re\{\lambda_1\}$','Interpreter','Latex');
grid on;

figure(11);
plot(q, real(lambda2));
title('$Re\{\lambda_2$\} vs $q$ for Case 3','Interpreter','Latex');
xlabel('$q$','Interpreter','Latex');
ylabel('$Re\{\lambda_2\}$','Interpreter','Latex');
grid on;

omega1 = sqrt(lambda1);
omega2 = sqrt(lambda2);

figure(12);
hold on
plot(q, abs(imag(omega1)));
plot(q, abs(imag(omega2)));
title('$\omega_1$ and $\omega_2$ vs $q$ for Case 3','Interpreter','Latex');
xlabel('$q$','Interpreter','Latex');
ylabel('$\omega$','Interpreter','Latex');
legend('$\omega_1$', '$\omega_2$', 'Interpreter', 'Latex');
grid on;

[maximum, index] = max(abs(imag(omega1)));
q_f = q(index)
U_f = sqrt(2*q_f/rho)

% Case 4
fprintf("Case 4: \n");
Sa = -1.5;

% a
q_d = k2/(e*S*dCLda)
U_d = sqrt(2*q_d/rho)

% b
q = 0:0.1:1000;

lambda1 = zeros(length(q), 1);
lambda2 = lambda1;
M = [m -Sa; -Sa Ia];

for i=1:length(q)
    
    K_eff = [k1 -q(i)*S*dCLda; 0 k2 - e*q(i)*S*dCLda];

    [X, D] = eig(K_eff, M);
    lambda1(i) = D(1, 1);
    lambda2(i) = D(2, 2);

    lambda1(i) = sqrt(lambda1(i)*-1);
    lambda2(i) = sqrt(lambda2(i)*-1);
end

figure(13);
plot(q, real(lambda1));
title('$Re\{\lambda_1$\} vs $q$ for Case 4','Interpreter','Latex');
xlabel('$q$','Interpreter','Latex');
ylabel('$Re\{\lambda_1\}$','Interpreter','Latex');
grid on;

figure(14);
plot(q, real(lambda2));
title('$Re\{\lambda_2$\} vs $q$ for Case 4','Interpreter','Latex');
xlabel('$q$','Interpreter','Latex');
ylabel('$Re\{\lambda_2\}$','Interpreter','Latex');
grid on;

omega1 = sqrt(lambda1);
omega2 = sqrt(lambda2);

figure(15);
hold on
plot(q, abs(imag(omega1)));
plot(q, abs(imag(omega2)));
title('$\omega_1$ and $\omega_2$ vs $q$ for Case 4','Interpreter','Latex');
xlabel('$q$','Interpreter','Latex');
ylabel('$\omega$','Interpreter','Latex');
legend('$\omega_1$', '$\omega_2$', 'Interpreter', 'Latex');
grid on;

[maximum, index] = max(abs(imag(omega1)));
q_f = q(index)
U_f = sqrt(2*q_f/rho)

%% Problem 3

x_init = [0 0 0 0.1];
tspan = [0 5];

[t, X] = ode45(@ode, tspan, x_init);

q1 = X(:,1);
q2 = X(:,2);

figure(16);
plot(t, q1);
title("Plot of q_1");
xlabel("Time (s)");
ylabel("Displacement");
title('$q_1$ vs $t$ for $q$ = 1.03$q_F$','Interpreter','Latex');
xlabel('$t$','Interpreter','Latex');
ylabel('$q_1$','Interpreter','Latex');
grid on

figure(2);
plot(t, q2);
title("Plot of q_2");
xlabel("Time (s)");
ylabel("Displacement");
title('$q_2$ vs $t$ for $q$ = 1.03$q_F$','Interpreter','Latex');
xlabel('$t$','Interpreter','Latex');
ylabel('$q_2$','Interpreter','Latex');
grid on

function dxdt = ode(t, x)
    m = 8.7;
    Sa = 5.4;
    Ia = 84;
    k1 = 6.3e4;
    k2 = 1.3e6;
    c = 5;
    e = 0.19;
    S = 150;
    dCLda = 1.87*pi;
    rho = 0.002377;
    q = 1.03*275.6;

    syms Q1 Q2
    L = q*S*dCLda*Q2;
    Mea = e*L;
    f = 1.3e6*(1+4*(Q2^2))*Q2;

    M = [m -Sa; -Sa Ia];
    F = [L-k1*Q1; Mea - f];

    qddot = M\F;

    % x = [q1 q2 q1dot q2dot]

    q1 = x(1);
    q2 = x(2);
    q1dot = x(3);
    q2dot = x(4);

    dxdt = zeros(size(x));
    dxdt(1) = q1dot;
    dxdt(2) = q2dot;
    dxdt(3) = subs(qddot(1), [Q1, Q2], [q1, q2]);
    dxdt(4) = subs(qddot(2), [Q1, Q2], [q1, q2]);
end
