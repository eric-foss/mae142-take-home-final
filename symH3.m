syms r11 r12 r13 r21 r22 r23 r31 r32 r33 r41 r42 r43

X_s1 = [8000; 100; 800; 0; 0; 0];
X_s2 = [100; 8000; 800; 0; 0; 0];
X_s3 = [1000; 2000; 5000; 0; 0; 0];

% X_s1 = [r11; r12; r13; 0; 0; 0];
% X_s2 = [r21; r22; r23; 0; 0; 0];
% X_s3 = [r31; r32; r33; 0; 0; 0];


syms x y z u v w t

X_sym = [x; y; z; u; v; w; t];

X0 = [5; 0; 0; 0; 0.26; 0; 0];



rho_s1 = sqrt((X_sym(1)-X_s1(1))^2 + (X_sym(2)-X_s1(2))^2 + (X_sym(3)-X_s1(3))^2) + t;
rho_s2 = sqrt((X_sym(1)-X_s2(1))^2 + (X_sym(2)-X_s2(2))^2 + (X_sym(3)-X_s2(3))^2) + t;
rho_s3 = sqrt((X_sym(1)-X_s3(1))^2 + (X_sym(2)-X_s3(2))^2 + (X_sym(3)-X_s3(3))^2) + t;

G = [rho_s1; rho_s2; rho_s3];

matlabFunction(G, 'File', 'computeG', "Vars", [x y z t]);
matlabFunction(jacobian(G, X_sym), 'File', 'computeH', "Vars", [x y z]);
