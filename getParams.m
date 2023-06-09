clear;
%Get constants from PDM integration
cd ./PDM_integration/
raw = readmatrix("/parameters.xlsx");
cd ../
MOI = raw(1:3, 2:4);
l = raw(4, 2);
m = raw(5, 2);
g = 9.81;

cd ./helper_files/

%Enumerate operating points for trajectory below:
q1 = angle2quat(0, 0, 0, 'ZYX')';
q1 = q1(2:end);
x1 = [10; 10; 0; 0; 0; 0; q1; 0; 0; 0];
u1 = [0; 0; m*g; 0];

q2 = angle2quat(0, 0, 0, 'ZYX')';
q2 = q2(2:end);
x2 = [5; 5; 0; 0; 0; 0; q2; 0; 0; 0];
u2 = [0; 0; m*g; 0];

q3 = angle2quat(0, 0, 0, 'ZYX')';
q3 = q3(2:end);
x3 = [1; 0; 0; 0; 0; 0; q3; 0; 0; 0];
u3 = [0; 0; m*g; 0];

inputLimits = [-pi/12, pi/12; -pi/12, pi/12; 0, 100; -2, 2];
throttleConsts = [-0.000112; 0.02072; -.268];

betaInputDelay = 0.001;
gammaInputDelay = 0.001;

%SIMULATE and create trajectory
fprintf("Creating Trajectory\n");
[x_set, u_set, t_set_unfiltered, Kset, tSegs] = get_trajectory(10000, 10, [x1, x2, x3], [u1, u2, u3], [m; l; g], MOI, inputLimits, throttleConsts);
[t_set, ia, ~] = unique(t_set_unfiltered, 'stable');
x_set = x_set(:, ia);
u_set = u_set(:, ia);
plotTrajectory(x_set, u_set, t_set, 0.75, 500, [m; l; g]);
startTime = t_set(1);
stopTime = t_set(end);

fprintf("Done initializing!\n");
cd ../