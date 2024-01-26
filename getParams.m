clear;
%Get constants from PDM integration
addpath('lib');
addpath('PDM_integration');
raw = readmatrix("/parameters.xlsx");
MOI = raw(1:3, 2:4);
l = raw(4, 2);
m = raw(5, 2);
g = 9.81;

%Enumerate operating points for trajectory below:
initEul = [pi/3, pi/12, pi/12];
q1 = angle2quat(initEul(1), initEul(2), initEul(3), 'XYZ')';
q1 = q1(2:end);
x1 = [10; 5; 5; 0; 0; 0; q1; 0; 0; 0];
u1 = [0; 0; m*g; 0];

q2 = angle2quat(0, 0, 0, 'XYZ')';
q2 = q2(2:end);
x2 = [5; 0; 0; 0; 0; 0; q2; 0; 0; 0];
u2 = [0; 0; m*g; 0];

q3 = angle2quat(0, 0, 0, 'XYZ')';
q3 = q3(2:end);
x3 = [0; 0; 0; 0; 0; 0; q3; 0; 0; 0];
u3 = [0; 0; m*g; 0];

%Maximum and minimum values for input for simulation
inputLimits = [-pi/12, pi/12; -pi/12, pi/12; 0, 100; -2, 2];

%Throttle -> Thrust force equation constants
throttleConsts = [-0.000112; 0.02072; -.268];

%Delay for servo angle inputs
betaInputDelay = 0.001;
gammaInputDelay = 0.001;

x_crit = [x1, x2, x3];
u_crit = [u1, u2, u3];

discretize_interval = 0.0025; % loop time in seconds

%SIMULATE and create trajectory
fprintf("Creating Trajectory\n");
[x_set, u_set, t_set, Kset, tSegs, startTime, stopTime] = get_trajectory(10000, 15, x_crit, u_crit, [m; l; g], MOI, inputLimits, throttleConsts);

%plot trajectory
plotTrajectory(x_set, u_set, t_set, 0.75, 300, [m; l; g]);

fprintf("Done initializing!\n");