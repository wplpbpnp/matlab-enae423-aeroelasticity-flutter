% Benjamin Stutzke
% ENAE 423
% Homework 8

%% Problem 3
% Free vibration of slender body
syms s real

E=1;
Area=1;
L=1;
m=1;
alpha = 1.5;

nedof = 2;
net = 4;
ntdof = net+1;
iconm = [1 2;2 3;3 4;4 5];

elength = L/net;
factorm = m*elength/6;
factork = E*Area/elength;
em = [2 1; 1 2]*factorm;
ek = [1 -1; -1 1]*factork;

gm = zeros(ntdof, ntdof);
gk = zeros(ntdof, ntdof);
gf = zeros(ntdof, 1);

nrdof = 5;
nreduced = [1 2 3 4 5];
gmr = zeros(nrdof, nrdof);
gkr = zeros(nrdof, nrdof);
gfr = zeros(nrdof, 1);

for lnum = 1:net
    iconv(:) = iconm(lnum, :);
    gk(iconv(:), iconv(:)) = gk(iconv(:), iconv(:))+ek(:,:);
    gm(iconv(:), iconv(:)) = gm(iconv(:), iconv(:))+em(:,:);
end

gf(1) = 1;
gk(:, 1) = gk(:, 1) + (alpha*E*Area/L)*gf;

fprintf("Global Stiffness Matrix:\n");
disp(gk);
fprintf("Global Mass Matrix:\n");
disp(gm);

gkr(:,:) = gk(nreduced(:), nreduced(:));
gmr(:,:) = gm(nreduced(:), nreduced(:));

[vec, p] = eig(gkr, gmr);
gv=zeros(ntdof,nrdof);
gv(nreduced(:),:)=vec(:,:);

nmode = 5;
msign = [1 -1 1 1 1];

beta = [0.9882, 3.5422, 6.5097, 9.5801, 12.684];

for i=1:nmode
    fprintf("\nMode %d\n", i);
    fprintf('\nExact Natural Frequency = %.4f sqrt(EA/mL^2)\n', beta(i));
    fprintf('\nNatural Frequency = %.4f sqrt(EA/mL^2)\n', vpa(sqrt(p(i,i))));

    figure;
    hold on

    N(1) = 1-s;
    N(2) = s;
    N = [N(1) N(2)];

    np = 500;

    for j=1:np
        s1 = (j-1)/np;
        N1 = subs(N, s, s1);

        u1(j) = (N1(1)*gv(1,i)+N1(2)*gv(2,i));
        u2(j) = (N1(1)*gv(2,i)+N1(2)*gv(3,i));
        u3(j) = (N1(1)*gv(3,i)+N1(2)*gv(4,i));
        u4(j) = (N1(1)*gv(4,i)+N1(2)*gv(5,i));
    end

    uall = [u1 u2 u3 u4];
    uabs = abs(uall);
    umax = max(uabs);
    uscaled = uall/umax;
    uscaled = msign(i) * uscaled;

    nptotal = np*net;
    delx = 1/nptotal;
    endpoint = 1-delx;
    xbar = 0:delx:endpoint;

    uexact = cos(beta(i)*(1-xbar/L));

    plot(xbar, uexact, 'r', 'LineWidth', 1);
    plot(xbar, uscaled, '--k', 'LineWidth', 1);
    legend('Exact', 'FEA');
    xlabel("x/L");
    title(sprintf("Mode %d", i));
    grid on
    hold off
end