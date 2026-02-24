% Benjamin Stutzke
% Midterm 1

%% Probem 1
K = [6 -2; -2 4];
M = [4 0; 0 2];

[X, D] = eig(K, M);

% D should contain all omega^2 values on diag
% X contains eigenvectors as columns

p_bar = diag(D);
freq = sqrt(p_bar) % natural frequency = freq*sqrt(k/m)

X_norm = bsxfun(@rdivide, X, max(abs(X)))

% Part C
phi_1 = X_norm(:, 1);
phi_2 = X_norm(:, 2);

phi_1'*K*phi_2
phi_1'*M*phi_2