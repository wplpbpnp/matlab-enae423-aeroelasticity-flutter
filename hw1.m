% Benjamin Stutzke
% ENAE 423
% Homework 1

%% Problem 1
tspan = 0:0.01:3;
c2 = 5.079;
omega_n = 4.655;
y = c2.*tspan.*exp(-omega_n.*tspan);
figure(1);
plot(tspan, y);
title("Problem 1c: Recoil displacement vs. time")
ylabel("Displacement (m)");
xlabel("Time (s)");

%% Problem 2
% Part A
global m;
global k;
global F0;
global OMEGA;
global c;

m = 1;          % kg
k = 25*(pi^2);  % N/m
F0 = 0.1;       % N
omega_n = sqrt(k/m);
c = 0;

% For 0.75
OMEGA = 0.75*omega_n;

tspan = [0 60];
xinit = [0 0];

[t, x] = ode45(@ode, tspan, xinit);

y = x(:,1);
ydot = x(:,2);

figure(2);
plot(t, y);
title("Problem 2a: Solutions for no damping at 0.75*\omega_n");
xlabel("Time (s)");
ylabel("Displacement (m)");

% For 0.95
OMEGA = 0.95*omega_n;
[t, x] = ode45(@ode, tspan, xinit);
y = x(:,1);
ydot = x(:,2);

figure(3);
plot(t, y);
title("Problem 2a: Solutions for no damping at 0.95*\omega_n");
xlabel("Time (s)");
ylabel("Displacement (m)");

% For 1.03
OMEGA = 1.03*omega_n;
[t, x] = ode45(@ode, tspan, xinit);
y = x(:,1);
ydot = x(:,2);

figure(4);
plot(t, y);
title("Problem 2a: Solutions for no damping at 1.03*\omega_n");
xlabel("Time (s)");
ylabel("Displacement (m)");

% Part B
zeta = 0.001;
c = 2*zeta*omega_n*m;

% For 0.75
OMEGA = 0.75*omega_n;

tspan = [0 60];
xinit = [0 0];

[t, x] = ode45(@ode, tspan, xinit);

y = x(:,1);
ydot = x(:,2);

figure(5);
plot(t, y);
title("Problem 2b: Solutions for damping ratio \zeta = 0.001 at 0.75*\omega_n");
xlabel("Time (s)");
ylabel("Displacement (m)");

% For 0.95
OMEGA = 0.95*omega_n;
[t, x] = ode45(@ode, tspan, xinit);
y = x(:,1);
ydot = x(:,2);

figure(6);
plot(t, y);
title("Problem 2b: Solutions for damping ratio \zeta = 0.001 at 0.95*\omega_n");
xlabel("Time (s)");
ylabel("Displacement (m)");

% For 1.03
OMEGA = 1.03*omega_n;
[t, x] = ode45(@ode, tspan, xinit);
y = x(:,1);
ydot = x(:,2);

figure(7);
plot(t, y);
title("Problem 2b: Solutions for damping ratio \zeta = 0.001 at 1.03*\omega_n");
xlabel("Time (s)");
ylabel("Displacement (m)");

function dxdt = ode(t, x)
    global m;
    global F0;
    global OMEGA;
    global c;
    global k;

    y = x(1);
    v = x(2);

    dxdt = zeros(size(x));
    dxdt(1) = v;
    dxdt(2) = (1/m)*(F0*sin(OMEGA*t)-(c*v+k*y));
end


