%% MAE142 Take-home Final Problem 1
% Eric Foss
clear; close all; clc;

dt = 1; %discrete time step [sec]

R = 10; %km
H = 5; %km
omega = 0.026; %rotation rate [rad/s]

s_bar = 0.26; %km/s

X0 = [R; 0; 0; 0; s_bar; 0];
C0 = [0.001*eye(3) zeros(3); zeros(3) 0.0001*eye(3)];

sigma_hat = 0.005;
sigma = sqrt(dt)*sigma_hat;
W = [zeros(3); sigma*eye(3)]; %noise matrix
C_w = eye(3); %covariance of noise


%% Part A B C

% Find A
A_hat = [zeros(3) eye(3); -omega^2*eye(3) zeros(3)]; %A_hat(6, 3) = 0;
A = eye(6) + dt*A_hat;

% propogate to find nominal
[t, X_cont] = ode45(@dynamics, 0:2*pi/omega, X0);

% Initialize data matrices
X_model = zeros(floor(2*pi/omega/dt)+1, 6); %model (no noise)
X_model(1, :) = X0';

X_rand = zeros(floor(2*pi/omega/dt)+1, 6); %with noise (true, i guess))
X_rand(1, :) = X0' + (sqrt(C0)*rand(6, 1))';

C = zeros(6, 6, floor(2*pi/omega/dt)+1); %covariance
C(:, :, 1) = C0;


for i = 2:((2*pi/omega/dt)+1)
    X_model(i, :) = (A*X_model(i-1, :)')';
    X_rand(i, :) = (A*X_rand(i-1, :)' + W*randn(3, 1))';
    C(:, :, i) = A*C(:, :, i-1)*A' + W*C_w*W';

end

figure(1); hold on;
plot(X_cont(:, 1), X_cont(:, 2), 'b-');
plot(X_model(:, 1), X_model(:, 2), 'ro');
xlabel('r_1 [km]'); ylabel('r_2 [km]'); 
legend('Continuous', sprintf('Discretized, dt = %0.1f', dt));
axis equal;

figure(2); hold on;
plot(X_cont(:, 1), X_cont(:, 2), 'b-');
plot(X_model(:, 1), X_model(:, 2), 'ro');
plot(X_rand(:, 1), X_rand(:, 2), 'gx');
xlabel('r_1 [km]'); ylabel('r_2 [km]'); 
legend('Continuous', sprintf('Model, dt = %0.1f', dt), sprintf('True, dt = %0.1f', dt));
axis equal

%% Part D

X_mc = zeros(floor(2*pi/omega/dt)+1, 6, 1000);

for k = 1:1000

    X_mc(1, :, k) = X0';

    for i = 2:((2*pi/omega/dt)+1)

        X_mc(i, :, k) = (A*X_mc(i-1, :, k)' + W*randn(3, 1))';

    end

end

figure(3);
subplot(1, 2, 1); hold on;
plot3(X_mc(:, 1, 25), X_mc(:, 2, 25), X_mc(:, 3, 25), 'bo');
plot3(X_mc(:, 1, 100), X_mc(:, 2, 100), X_mc(:, 3, 100), 'r.');
xlabel('r_1 [km]'); ylabel('r_2 [km]'); zlabel('r_3 [km]');
legend('X_{25}', 'X_{100}');
view(30, 15);

subplot(1, 2, 2); hold on;
plot3(X_mc(:, 4, 25), X_mc(:, 5, 25), X_mc(:, 6, 25), 'bo');
plot3(X_mc(:, 4, 100), X_mc(:, 5, 100), X_mc(:, 6, 100), 'r.');
xlabel('v_1 [km/s]'); ylabel('v_2 [km/s]'); zlabel('v_3 [km/s]');
legend('X_{25}', 'X_{100}');
view(30, 15);




