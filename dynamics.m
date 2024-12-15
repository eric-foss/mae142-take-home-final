function xdot = dynamics(~, x, const)

omega = 0.026;

xdot = [x(4:6); -omega^2*x(1:2); 0];