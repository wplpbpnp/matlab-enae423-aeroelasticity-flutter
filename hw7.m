% Benjamin Stutzke
% ENAE 423
% Homework 7

%% Problem 1
m = 1;
L = 1;
EI0 = 1;

syms s

N_array = {[s^2 s^3]; [s^2 s^3 s^4]; [s^2 s^3 s^4 s^5]};

for prob=1:length(N_array)
    N = N_array{prob};

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
    end
end
%% Problem 2
m = 1;
L = 1;
EIy = 1;

syms s

% check assumptions by assuming different N's
N_array = {[1 s s^2 s^3 s^4 s^5];}; %[1 s s^2 s^3 s^4 s^5 s^6]; [1 s s^2 s^3]; [1 s s^2 s^3 s^4]; [1 s s^2 s^3 s^4 s^5 s^6 s^7 s^8 s^9]};
sign_array = {};

for index = 1:length(N_array)
    sign_array{index} = ones(size(N_array{index}));
end

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

        np = 50;

        for j=1:np
            s1 = j/np;
            newN = subs(N, s, s1);
            w1(j) = dot(newN, X(:, k));
        end

        msign = sign_array{index};
        wall=[0.0 w1]; % Unscaled mode
        wabs=abs(wall);
        wmax=max(wabs);
        wscaled=wall/wmax; % Scaled mode
        wscaled=msign(k)*wscaled; % Adjust the sign of mode # i
        delx=1/np;
        xbar=0:delx:1;
        figure(k);
        plot(xbar,wscaled,'--k','LineWidth',1)
        xlabel('x/L')
        title(sprintf('Mode %d',k))
        grid on
        hold off
    end
end