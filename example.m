% Chapter 7 Problem 7.11 Part (a)
% Assumed mode method
% Two-term approximation for cantilevered beam
clear; close all; clc
format short e
m = 1; % Mass per length
L = 1; % Total length
EI0 = 1; % Bending rigidity
syms s
N1 = s^2;
N2 = s^3;
N = [N1 N2];
factorM = m*L;
igrandM = N.'*N*(1-s/2);
M = int(igrandM,s,0,1)*factorM; % Mass matrix
M = double(M);
fprintf('Nondimensional Mass Matrix:\n')
disp(M)
factorK = EI0/(L^3);
B =diff(N,s,2);
igrandK = B.'*B*(1-s/2);
K = int(igrandK,s,0,1)*factorK; % Stiffness matrix
K = double(K);
fprintf('Nondimensional Stiffness Matrix:\n')
disp(K)
% K*x=Lambda*M*x where Lambda = omega^(mL^4)/(2*EIy)
[gv,p]=eig(K,M);
nmode=2; % The number of modes to be reported
msign = [-1 1]; % Mode sign control. Assign 1 intially for all modes.
% Check
% the signs of each mode, and then assign 1 or -1 to adjust the max value
% of each mode to be equal to 1.
for i=1:nmode
fr(i)=sqrt(p(i,i));
fprintf('\nMode %d\n',i)
fprintf('\nNatural Frequency = %.4f sqrt(EIy/mL^4)\n',fr(i))
figure
hold on
% Modes
N1 = s^2;
N2 = s^3;
N = [N1 N2];
np=50;
for j=1:np
s1=j/np;
NN=subs(N,s,s1);
w1(j)=(NN(1)*gv(1,i)+NN(2)*gv(2,i));
w2(j) = dot(NN, gv(:,i));
end

wall=[0.0 w1]; % Unscaled mode
wabs=abs(wall);
wmax=max(wabs);
wscaled=wall/wmax; % Scaled mode
wscaled=msign(i)*wscaled; % Adjust the sign of mode # i
delx=1/np;
xbar=0:delx:1;
plot(xbar,wscaled,'--k','LineWidth',1)
legend('Two mode solution')
xlabel('x/L')
title(sprintf('Mode %d',i))
grid on
hold off
end