cd ./PDM_integration/
raw = readmatrix("/parameters.xlsx");
cd ../
MOI = raw(1:3, 2:4);
l = raw(4, 2);
m = raw(5, 2);

g = 9.81;

q1 = angle2quat(0, 0, 0, 'ZYX')';
q1 = q1(2:end);
x1 = [0; 0; 0; 0; 0; 0; q1; 0; 0; 0];
u1 = [0; 0; m*g; 0];

q2 = angle2quat(0, 0, 0, 'ZYX')';
q2 = q2(2:end);
x2 = [10; 2.5; 0; 0; 0; 0; q2; 0; 0; 0];
u2 = [0; 0; m*g; 0];

q3 = angle2quat(0, 0, 0, 'ZYX')';
q3 = q3(2:end);
x3 = [0; 5; 0; 0; 0; 0; q3; 0; 0; 0];
u3 = [0; 0; m*g; 0];

q4 = angle2quat(0, 0, 0, 'ZYX')';
q4 = q4(2:end);
x4 = [1; 10; 0; 0; 0; 0; q4; 0; 0; 0];
u4 = [0; 0; m*g; 0];

q5 = angle2quat(0, 0, 0, 'ZYX')';
q5 = q5(2:end);
x5 = [5; 4.5; 0; 0; 0; 0; q5; 0; 0; 0];
u5 = [0; 0; m*g; 0];

q6 = angle2quat(0, 0, 0, 'ZYX')';
q6 = q6(2:end);
x6 = [10; 10; 0; 0; 0; 0; q6; 0; 0; 0];
u6 = [0; 0; m*g; 0];

rmax = [1; 1; 1];
vmax = [1; 1; 1];
qmax = angle2quat(pi/12, pi/12, pi/12, 'ZYX')';
qmax = qmax(2:end);
omegamax = [pi/2; pi/2; pi/2];
xmax = [rmax; vmax; qmax; omegamax];

umax = [pi/12; pi/12; m*g*5; 1];

%Cost function Q and R matrices
Q = diag(xmax.^-2);
R = diag(umax.^-2);

fprintf("Creating Trajectory\n");
[xset, uset, tset] = get_trajectory(10000, [x1, x2, x3, x4, x5, x6], [u1, u2, u3, u4, u5, u6], [m; l; g], MOI);
plotTrajectory(xset, uset, tset, 0.5, 200, [m; l; g]);

fprintf("Done initializing!\n");