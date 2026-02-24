% Benjamin Stutzke
% ENAE 423
% Homework 7

%% Problem 1
m = 1;
L = 1;
EI0 = 1;

syms s

N_array = {[s^2 s^3]; [s^2 s^3 s^4]; [s^2 s^3 s^4 s^5]};
sign_array = {[-1 1]; [1 1 1]; [1 1 1 1]};

for prob=1:length(N_array)
    N = N_array{prob};
    msign = sign_array{prob};

    fprintf("For N = ");
    disp(N);

    igrandM = N.'*N*(1-s/2);
    M = int(igrandM, s, 0, 1) * m*L;
    Mbar = double(M);
    
    B = diff(N, s, 2);
    igrandK = B'.*B*(1-s/2);
    K = int(igrandK, s, 0 ,1)*EI0/(L^3);
    Kbar = double(K);
    
    % KX = lambda*MX, with lambda = omega^2*(mL^4)/(EIy)
    [X, D] = eig(Kbar, Mbar);
    
    nmode = length(N);
    
    for k=1:nmode
        freq(k) = sqrt(D(k,k));
        fprintf("Mode %d\n", k);
        fprintf("omega = %.4f sqrt(EIy/mL^4)\n", freq(k));
        fprintf("\n");

        figure
        hold on

        np = 500;
        np1 = np+1;

        for j=1:np1
            s1 = (j-1)/np;
            NN = subs(N, s, s1);
            w1(j) = NN(1, :)*X(:, k);
        end

        wall = w1;
        wabs = abs(wall);
        wmax = max(wabs);
        wscaled = wall/wmax;

        wscaled = msign(k)*wscaled;

        delx = 1/np;
        xbar = 0:delx:1;

        plot(xbar, wscaled, '--k', 'LineWidth', 1);
        xlabel('x/L');
        title(sprintf("%d Mode Solution: Mode %d", prob+1, k));
        grid on
        hold off
    end
end
clear
%% Problem 2
m = 1;
L = 1;
EIy = 1;

syms s

N_array = {[1 s s^2 s^3 s^4 s^5];};
sign_array = {[1 1 1 1 -1 1]};

for index=1:length(N_array)
    N = N_array{index};
    fprintf("For N = ");
    disp(N);    

    igrandM = N.'*N;
    M = int(igrandM, s, 0, 1) * m*L;
    Mbar = double(M);
    
    B = diff(N, s, 2);
    igrandK = B'.*B;
    K = int(igrandK, s, 0 ,1)*EIy/(L^3);
    Kbar = double(K);
    
     % KX = lambda*MX, with lambda = omega^2*(mL^4)/(EIy)
    [X, D] = eig(Kbar, Mbar);
    
    nmode = length(N);
    
    for k=1:nmode
        freq(k) = sqrt(D(k,k));
        fprintf("Mode %d\n", k);
        fprintf("omega = %.4f sqrt(EIy/mL^4)\n", freq(k));
        fprintf("\n");

        np = 500;

        for j=1:np+1
            s1 = (j-1)/np;
            newN = subs(N, s, s1);
            w1(j) = newN(1, :) * X(:, k);
        end

        msign = sign_array{index};
        wall=w1; % Unscaled mode
        wabs=abs(wall);
        wmax=max(wabs);
        wscaled=wall/wmax; % Scaled mode
        wscaled=msign(k)*wscaled; % Adjust the sign of mode # i
        delx=1/np;
        xbar=0:delx:1;
        figure(k);
        plot(xbar,wscaled,'--k','LineWidth',1)
        xlabel('x/L')
        title(sprintf('Mode %d', k))
        grid on
        hold off
    end
end