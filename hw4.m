% Benjamin Stutzke
% ENAE 423
% Homework 4

%% Problem 3
global BETA;
global K;
global M;

M = 1;          % kg
K = 4*(pi^2);   % N/m
A = 0.01;       % m

xinit = [A 0 0 0];

% beta = 1;
BETA = 1;
tspan = [0 10];

[t, x] = ode45(@ode, tspan, xinit);

q1 = x(:,1);
q2 = x(:,3);

figure(1);
plot(t, q1./A);
title("Plot of q_1/A for \beta = 1.0");
xlabel("Time (s)");
ylabel("Displacement");

figure(2);
plot(t, q2./A);
title("Plot of q_2/A for \beta = 1.0");
xlabel("Time (s)");
ylabel("Displacement");

% beta = 0.05;
BETA = 0.05;
tspan = [0 100];

[t, x] = ode45(@ode, tspan, xinit);

q1 = x(:,1);
q2 = x(:,3);

figure(3);
plot(t, q1./A);
title("Plot of q_1/A for \beta = 0.05");
xlabel("Time (s)");
ylabel("Displacement");

figure(4);
plot(t, q2./A);
title("Plot of q_2/A for \beta = 0.05");
xlabel("Time (s)");
ylabel("Displacement");

function dxdt = ode(t,x)

    global BETA;
    global K;
    global M;

    % X = [q1, q1dot, q2, q2dot]
    q1 = x(1);
    q1dot = x(2);
    q2 = x(3);
    q2dot = x(4);

    dxdt = zeros(size(x));
    dxdt(1) = q1dot;
    dxdt(2) = (BETA*K*q2 - (K+BETA*K)*q1)/M;
    dxdt(3) = q2dot;
    dxdt(4) = (BETA*K*q1 - (K+BETA*K)*q2)/M;
end

