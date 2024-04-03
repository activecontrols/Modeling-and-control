clear;
close all;
%Get constants from PDM integration
addpath('lib');
addpath('PDM_integration');
raw = readmatrix("/parameters.xlsx");
MOI = raw(1:3, 2:4);
l = raw(4, 2);
m = raw(5, 2);
g = 9.81;

%Enumerate operating points for trajectory below:
% initEul = [pi/3, pi/12, pi/12];
initEul = [0, 0, 0];
% order is roll, pitch, yaw?
q1 = angle2quat(initEul(1), initEul(2), initEul(3), 'XYZ')'; % todo: sim depends on q1 rn
q1 = q1(2:end);
x1 = [10; 0; 0; 0; 0; 0; q1; 0; 0; 0]; % todo: sim depends on x1 rn
u1 = [0; 0; m*g; 0];

q2 = angle2quat(0, 0, 0, 'XYZ')';
q2 = q2(2:end);
x2 = [15; 15; -15; 0; 0; 0; q2; 0; 0; 0];
u2 = [0; 0; m*g; 0];

q3 = angle2quat(0, 0, 0, 'XYZ')';
q3 = q3(2:end);
x3 = [1; 0; 0; 0; 0; 0; q3; 0; 0; 0];
u3 = [0; 0; m*g; 0];

%Maximum and minimum values for input for simulation
inputLimits = [-pi/12, pi/12; -pi/12, pi/12; 0, 100; -2, 2];

%Throttle -> Thrust force equation constants
% throttleConsts = [-0.000112; 0.02072; -.268];
throttleConsts = [71.0072; 0.0504; 67.9056];

% throttleA = 71.0072;
% throttleB = 0.0504;
% throttleC = 67.9056;

%Delay for servo angle inputs
betaInputDelay = 0.001;
gammaInputDelay = 0.001;

x_crit = [x1, x2, x3];
u_crit = [u1, u2, u3];

discretize_interval = 0.0025; % loop time in seconds

%SIMULATE and create trajectory
fprintf("Creating Trajectory\n");
[x_set, u_set, t_set, Kset, tSegs, startTime, stopTime] = get_trajectory(10000, 15, x_crit, u_crit, [m; l; g], MOI, inputLimits, throttleConsts);

% figure(1)

%plot trajectory
plotTrajectory(x_set, u_set, t_set, 0.75, 300, [m; l; g]);

% sets the model workspace
constructSimInput(MOI, l, m, g, inputLimits, throttleConsts, betaInputDelay, gammaInputDelay, discretize_interval, ...
    x_set, u_set, t_set, Kset, tSegs, startTime, stopTime, q1, x1);

fprintf("Set simulation params!\n");

fprintf("Done initializing!\n");

fprintf("Running simulation\n")

clear;

sim("Model_Mk2.slx")

xzy = yout(:, 5:7);

hold on

plot3(xzy(:, 1), xzy(:, 3), xzy(:, 2), 'LineWidth', 2) % plot trajectory from sim.

hold off