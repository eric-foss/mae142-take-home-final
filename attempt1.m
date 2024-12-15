%% MAE142 Take-home Final Problem 2
% Eric Foss
clear; close all; clc;

dt = 4;

R = 10;
omega = 0.026;

s_bar = 0.26;

X0 = [R/2; 0; 0; 0; s_bar; 0];
C0 = [10*eye(3) zeros(3, 4); zeros(4, 3) 1*eye(4)]; C0(7, 7) = 100;

sigma_hat = 0.005;
sigma = sqrt(dt)*sigma_hat;
W = [zeros(3); sigma*eye(3)];
C_w = eye(3);

sigma_nu = sqrt(1e-4);
C_nu = sigma_nu^2*eye(4);

X_s1 = [8000; 100; 800; 0; 0; 0];
X_s2 = [100; 8000; 800; 0; 0; 0];
X_s3 = [1000; 2000; 5000; 0; 0; 0];
X_s4 = [2000; 100; 2000; 0; 0; 0];

%% Configuration    


A_hat = [zeros(3) eye(3); -omega^2*eye(3) zeros(3)];
A = eye(6) + dt*A_hat;
 

X_real = zeros(floor(2*pi/omega/dt)+1, 6);
X_real(1, :) = X0';

X_model = zeros(floor(2*pi/omega/dt)+1, 7);
X_model(1, :) = [X0' 0];

C_model = zeros(6, 6, floor(2*pi/omega/dt)+1);
C(:, :, 1) = C0(1:6, 1:6);

X_bar = zeros(floor(2*pi/omega/dt)+1, 7);
X_bar(1, :) = [X0' 0];

X_est = zeros(floor(2*pi/omega/dt)+1, 7);
X_est(1, :) = [X0' 0];

M = zeros(7, 7, floor(2*pi/omega/dt)+1);
M(:, :, 1) = C0;

P = zeros(7, 7, floor(2*pi/omega/dt)+1);
P(:, :, 1) = C0;

H = zeros(4, 7, floor(2*pi/omega/dt)+1);
H(:, :, 1) = computeH(X0(1), X0(2), X0(3));

obsv = zeros(floor(2*pi/omega/dt)+1, 4);
obsv(1, :) = (computeG(X0(1), X0(2), X0(3), 0) + sigma_nu*randn(4, 1));

for i = 2:((2*pi/omega/dt)+1)

    X_real(i, :) = (A*X_real(i-1, :)' + W*randn(3, 1))';

    X_model(i, 1:6) = (A*X_model(i-1, 1:6)')'; X_model(i, 7) = X_model(i-1, 7)+dt;
    C_model(:, :, i)= A*C_model(:, :, i-1)*A' + W*C_w*W';

    % Dynamics Update
    X_bar(i, 1:6) = (A*X_est(i-1, 1:6)')'; X_bar(i, 7) = X_est(i-1, 7);
    M(1:6, 1:6, i) = A*P(1:6, 1:6, i-1)*A' + W*C_w*W'; M(7, 7, i) = 100;

    % Observations
    obsv(i, :) = (computeG(X_real(i, 1), X_real(i, 2), X_real(i, 3), dt*(i-1)) + sigma_nu*randn(4, 1))';

    H(:, :, i) = computeH(X_model(i, 1), X_model(i, 2), X_model(i, 3));

    %P(:, :, i) = inv(inv(M(:, :, i)) + H(:, :, i)'*inv(C_nu)*H(:, :, i));

    K = M(:, :, i)*H(:, :, i)'*inv(H(:, :, i)*M(:, :, i)*H(:, :, i)' + C_nu);

    %X_est(i, :) = (X_model(i, :)' + P(:, :, i)*H(:, :, i)'*inv(C_nu)*(obsv(i, :)' - H(:, :, i)*X_model(i, :)'))';

    X_est(i, :) = (X_model(i, :)' + K*(obsv(i, :)' - computeG(X_model(i, 1), X_model(i, 2), X_model(i, 3), (i-1)*dt) - H(:, :, i)*X_model(i, :)'))';
    P(:, :, i) = (eye(7) - K*H(:, :, i))*M(:, :, i);

    %computeG(X_model(i, 1), X_model(i, 2), X_model(i, 3), X_model(i, 7)) 
end


%% Figures

figure(1); hold on;
plot3(X_real(:, 1), X_real(:, 2), X_real(:, 3), 'ob');
plot3(X_est(:, 1), X_est(:, 2), X_est(:, 3), 'r.');

