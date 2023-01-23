function [traj, t] = calc_trajectory(xInitial, xFinal, xWaypoints, numPoints, endTime)

posInit = xInitial(1:3);
posFinal = xFinal(1:3);
if ~isempty(xWaypoints)
    posWaypoints = xWaypoints(1:3, :);
else 
    posWaypoints = [];
end

posTraj = [posInit, posWaypoints, posFinal];

posTrajNew = interparc(numPoints, posTraj(1, :), posTraj(2, :), posTraj(3, :), 'spline')';
X1 = posTrajNew(1, :);
X2 = posTrajNew(2, :);
X3 = posTrajNew(3, :);
t = linspace(0, endTime, numPoints);
VX1 = gradient(X1, t);
VX2 = gradient(X2, t);
VX3 = gradient(X3, t);

subplot(2, 2, 1);
plot3(posTraj(1, :), posTraj(2, :), posTraj(3, :), 'r.', X1, X2, X3, 'b-');
subplot(2, 2, 2);
plot3(VX1, VX2, VX3, 'g-');

oriInit = xInitial(4:6);
oriFinal = xFinal(4:6);
if ~isempty(xWaypoints)
    oriWaypoints = xWaypoints(1:3, :);
else 
    oriWaypoints = [];
end

oriTraj = [oriInit, oriWaypoints, oriFinal];

oriTrajNew = interparc(numPoints, oriTraj(1, :), oriTraj(2, :), oriTraj(3, :), 'spline')';
phi = oriTrajNew(1, :);
theta = oriTrajNew(2, :);
psi = oriTrajNew(3, :);
omega_phi = gradient(phi, t);
omega_theta = gradient(theta, t);
omega_psi = gradient(psi, t);

subplot(2, 2, 3);
plot3(oriTraj(1, :), oriTraj(2, :), oriTraj(3, :), 'r.', phi, theta, psi, 'b-');
subplot(2, 2, 4);
plot3(omega_phi, omega_theta, omega_psi, 'g-');

traj = [X1; X2; X3; phi; theta; psi; VX1; VX2; VX3; omega_phi; omega_theta; omega_psi];
end

