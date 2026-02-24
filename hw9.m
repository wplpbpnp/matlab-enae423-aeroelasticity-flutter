% Benjamin Stutzke
% ENAE 423
% Homework 9

M = [1 0; 0 1];
C = [0.8 -0.8; -0.8 1.6];
K = [400 -400; -400 800];
Fapp = [0; 0];

tstep = 0.015;
tspan = 0:tstep:10;
InitCon = [0; 0; 1; 0];

[ts, xv] = ode45(@TimeRate, tspan, InitCon);

SNR = 150;
LL = length(xv(:,1));
noise1 = std(xv(:,1))/SNR*randn(LL,1);
noise2 = std(xv(:,2))/SNR*randn(LL,1);
q(:,1) = xv(:, 1) + noise1;
q(:,2) = xv(:, 2) + noise2;

figure(1);
hold on
plot(ts, q(:, 1))
plot(ts, q(:, 2))
title("q1 and q2 vs Time");
xlabel("Time (s)");
ylabel("Displacement (m)");

D01 = (q([1:length(q) - 4], 1))';
D02 = (q([1:length(q) - 4], 2))';

D0 = [D01; D02];

D1 = [(q([2:length(q) - 3],1))';(q([2:length(q) - 3],2))'];
D2 = [(q([3:length(q) - 2],1))';(q([3:length(q) - 2],2))'];
D3 = [(q([4:length(q) - 1],1))';(q([4:length(q) - 1],2))'];
D4 = [(q([5:length(q) - 0],1))';(q([5:length(q) - 0],2))'];

Dstar = [D0 D1; D1 D2];
Dstarhat = [D2 D3; D3 D4];

DstarT = Dstar';
inverse = inv(Dstar*DstarT);

A1 = Dstarhat*DstarT*inverse;
I = eye(4);

[x, D] = eig(A1, I);

lambda1 = log(D(1,1))/(2*tstep)
lambda2 = log(D(2,2))/(2*tstep)
lambda3 = log(D(3,3))/(2*tstep)
lambda4 = log(D(4,4))/(2*tstep)

freq1 = abs(imag(lambda1))
freq2 = abs(imag(lambda3))

function output = TimeRate(t, XV)

    M = [1 0; 0 1];
    C = [0.8 -0.8; -0.8 1.6];
    K = [400 -400; -400 800];
    Fapp = [0; 0];

    q1 = XV(1);
    q2 = XV(2);
    v1 = XV(3);
    v2 = XV(4);

    q = [q1; q2];
    qdot = [v1; v2];
    vdot = M\ (Fapp-K*q - C*qdot);

    output = [qdot(1); qdot(2); vdot(1); vdot(2)];
end