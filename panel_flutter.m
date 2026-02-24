% Benjamin Stutzke
% ENAE 423
% Final

m = 1; % Mass per length
L = 1; % Total length
EIy = 1; % Bending rigidity
syms s

nterm = 5;

N1=sin(pi*s);
N2=sin(2*pi*s);
N3=sin(3*pi*s);
N4=sin(4*pi*s);
N5=sin(5*pi*s);
N= [N1 N2 N3 N4 N5];
N_t=transpose(N);

factorM = m*L;
M = int(N_t*N,s,0,1)*factorM; % Mass matrix
M= double(M);
fprintf('Nondimensional Mass Matrix:\n')
disp(M)

factorK = EIy/(L^3);
B = diff(N,s,2);
B_t =transpose(B);
K = int(B_t*B,s,0,1)*factorK; % Bending Stiffness matrix
K=double(K);
fprintf('Nondimensional Stiffness Matrix:\n')
disp(K)

% Nondimesional nonsymmetric K_a matrix
A = diff(N,s);
Ka = -int(N_t*A,s,0,1);
Ka = double(Ka);
fprintf('Nondimensinal K_a Matrix:\n')
disp(Ka)

% Carry out eigenvalue analysis for given values of q_bar.
for j=1:451
    q_bar=1.0*(j-1);
    Keff=K-(q_bar*Ka);
    Lambda=eig(Keff,M);
    for k=1:nterm
        lamb(k)=sqrt(-Lambda(k));
        sigma(k)=real(lamb(k));
        omega(k)=abs(imag(lamb(k)));
    end
    omega =double(omega);
    [freq,II]=mink(omega,2);
    om1(j)=freq(1);
    om2(j)=freq(2);
    sig1(j)=sigma(II(1));
    sig2(j)=sigma(II(2));
end

xbar=0:1:450;

figure
hold on
plot(xbar,sig1,xbar,sig2,'r','LineWidth',1)
xlabel('$\bar{q}$','Interpreter','Latex','FontSize',14)
ylabel('$\bar{\sigma_k}$','Interpreter','Latex','FontSize',14)
grid on
hold off

figure
hold on
plot(xbar,om1,'--k','LineWidth',1)
plot(xbar,om2,'r','LineWidth',1)
legend('$\bar{\omega_1}$','$\bar{\omega_2}$','Interpreter','Latex','FontSize',14,'Location','NorthEast')
xlabel('$\bar{q}$','Interpreter','Latex','FontSize',14)
ylabel('$\bar{\omega_k}$','Interpreter','Latex','FontSize',14)
grid on
hold off
% According to the plots, critical q_bar is about 340.