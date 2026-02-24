% Benjamin Stutzke
% ENAE 423
% Homework 3

%% Problem 1
% Part A
k = [1 0 -1; 0 1 -1; -1 -1 3];
M = [1 0 0; 0 1 0; 0 0 1];

[X, D] = eig(k, M);

% D should contain all omega^2 values on diag
% X contains eigenvectors as columns

p_bar = diag(D);
freq = sqrt(p_bar) % natural frequency = freq*sqrt(k/m)

X_norm = bsxfun(@rdivide, X, max(abs(X)))

% Part C
phi_1 = X_norm(:, 1);
phi_2 = X_norm(:, 2);

phi_1'*k*phi_2
phi_1'*M*phi_2

%% Problem 2
% Part A
k = [1 0 -1; 0 1 -1; -1 -1 2];
M = [1 0 0; 0 1 0; 0 0 1];

[X, D] = eig(k, M);

% D should contain all omega^2 values on diag
% X contains eigenvectors as columns

p_bar = diag(D);
freq = sqrt(p_bar) % natural frequency = freq*sqrt(k/m)

X_norm = bsxfun(@rdivide, X, max(abs(X)))

% Part C
phi_1 = X_norm(:, 1);
phi_2 = X_norm(:, 2);

phi_1'*k*phi_2
phi_1'*M*phi_2