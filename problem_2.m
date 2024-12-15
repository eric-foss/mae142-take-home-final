%% MAE142 Take-home Final Problem 2
% Eric Foss
clear; close all;

num_sat = 4;

if num_sat == 4
    symH4;
else
    symH3;
end

dt = 4;

R = 10;
omega = 0.026;

s_bar = 0.26;

X0 = [R/2; 0; 0; 0; s_bar; 0];

simga_t = 10;
c10 = 10; c30 = 1;
C0 = [c10*eye(3) zeros(3, 4); zeros(4, 3) c30*eye(4)]; C0(7, 7) = simga_t^2;

A_hat = [zeros(3) eye(3); -omega^2*eye(3) zeros(3)];
A = eye(6) + dt*A_hat;

sigma_hat = 0.005;
sigma = sqrt(dt)*sigma_hat;
W = [zeros(3); sigma*eye(3)];
C_w = eye(3);

sigma_nu = sqrt(1e-4);
C_nu = sigma_nu^2*eye(num_sat);

X_s1 = [8000; 100; 800; 0; 0; 0];
X_s2 = [100; 8000; 800; 0; 0; 0];
X_s3 = [1000; 2000; 5000; 0; 0; 0];
X_s4 = [2000; 100; 2000; 0; 0; 0];

%% Configuration    

X_nom = zeros(floor(2*pi/omega/dt)+1, 7);
X_nom(1, :) = [X0' 0];

X = zeros(floor(2*pi/omega/dt)+1, 7);
X(1, :) = [X0' 0] + (sqrt(C0)*randn(7, 1))';

x_bar = zeros(floor(2*pi/omega/dt)+1, 7);
x_bar(1, :) = zeros(1, 7);

X_bar = zeros(floor(2*pi/omega/dt)+1, 7);
X_bar(1, :) = X_nom(1, :);

x_hat = zeros(floor(2*pi/omega/dt)+1, 7);
x_hat(1, :) = x_bar(1, :);

X_est = zeros(floor(2*pi/omega/dt)+1, 7);
X_est(1, :) = X_nom(1, :) + x_hat(1, :);

Y = zeros(floor(2*pi/omega/dt)+1, num_sat);
Y(1, :) = (computeG(X(1, 1), X(1, 2), X(1, 3), 0) + sigma_nu*rand(num_sat, 1));

y = zeros(floor(2*pi/omega/dt)+1, num_sat);
y(1, :) = (computeG(X(1, 1), X(1, 2), X(1, 3), 0) + sigma_nu*rand(num_sat, 1) - computeG(X_nom(1, 1), X_nom(1, 2), X_nom(1, 3), X_nom(1, 7)));

M = zeros(7, 7, floor(2*pi/omega/dt)+1);
M(:, :, 1) = C0;

H = zeros(num_sat, 7, floor(2*pi/omega/dt)+1);
H(:, :, 1) = computeH(X_nom(1, 1), X_nom(1, 2), X_nom(1, 3));

P = zeros(7, 7, floor(2*pi/omega/dt)+1);
P(:, :, 1) = inv(inv(M(:, :, 1)) + H(:, :, 1)'*inv(C_nu)*H(:, :, 1));

for i = 2:((2*pi/omega/dt)+1) 

    X(i, 1:6) = (A*X(i-1, 1:6)' + W*randn(3, 1))'; X(i, 7) = X(i-1, 7) + dt;

    X_nom(i, 1:6) = (A*X_nom(i-1, 1:6)')'; X_nom(i, 7) = X_nom(i-1, 7) + dt;

    X_bar(i, 1:6) = (A*X_est(i-1, 1:6)')'; X_bar(i, 7) = X_bar(i-1, 7) + dt;

    x_bar(i, 1:6) = (A*x_hat(i-1, 1:6)')'; x_bar(i, 7) = x_bar(i-1, 7) + dt;
    M(1:6, 1:6, i) = A*P(1:6, 1:6, i-1)*A' + W*C_w*W'; M(7, 7, i) = simga_t^2;

    Y(i, :) = (computeG(X(i, 1), X(i, 2), X(i, 3), X(i, 7)) + sigma_nu*randn(num_sat, 1))';

    y(i, :) = (Y(i, :)' - computeG(X_nom(i, 1), X_nom(i, 2), X_nom(i, 3), X_nom(i, 7)))';

    H(:, :, i) = computeH(X_nom(i, 1), X_nom(i, 2), X_nom(i, 3));

    P(:, :, i) = inv(inv(M(:, :, i)) + H(:, :, i)'*inv(C_nu)*H(:, :, i));

    x_hat(i, :) = (x_bar(i, :)' + P(:, :, i)*H(:, :, i)'*inv(C_nu)*(y(i, :)' - H(:, :, i)*x_bar(i, :)'))';

    X_est(i, :) = X_nom(i, :) + x_hat(i, :);

end


%% Figures

figure(1); hold on;
subplot(2, 2, 1);
plot(X(:, 7), Y(:, 1)-X(:, 7));
xlabel("Time [sec]"); ylabel("y^1 [sec]");
subplot(2, 2, 2);
plot(X(:, 7), Y(:, 2)-X(:, 7));
xlabel("Time [sec]"); ylabel("y^2 [sec]");
subplot(2, 2, 3);
plot(X(:, 7), Y(:, 3)-X(:, 7));
xlabel("Time [sec]"); ylabel("y^3 [sec]");
subplot(2, 2, 4);
if num_sat == 4
    plot(X(:, 7), Y(:, 4)-X(:, 7));
    xlabel("Time [sec]"); ylabel("y^4 [sec]");
end



figure(2); hold on;
plot3(X_nom(:, 1), X_nom(:, 2), X_nom(:, 3), 'bo');
plot3(X(:, 1), X(:, 2), X(:, 3), 'r.');
plot3(X_bar(:, 1), X_bar(:, 2), X_bar(:, 3), 'gx');
plot3(X_est(:, 1), X_est(:, 2), X_est(:, 3), 'kx');
legend("Model", "True", "A Priori", "A Posteriori");
xlabel("r_1 [km]"); ylabel("r_2 [km]"); zlabel("r_3 [km]");
view(30, 15);

figure(3); hold on;
plot3(X_nom(:, 4), X_nom(:, 5), X_nom(:, 3), 'bo');
plot3(X(:, 4), X(:, 5), X(:, 6), 'r.');
plot3(X_bar(:, 4), X_bar(:, 5), X_bar(:, 6), 'gx');
plot3(X_est(:, 4), X_est(:, 5), X_est(:, 6), 'kx');
legend("Model", "True", "A Priori", "A Posteriori");
xlabel("v_1 [km/s]"); ylabel("v_2 [km/s]"); zlabel("v_3 [km/s]");
view(30, 15);




