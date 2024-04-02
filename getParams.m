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
% initEul = [pi/3, pi/12, pi/12];
initEul = [0, 0, 0];
% order is roll, pitch, yaw?
q1 = angle2quat(initEul(1), initEul(2), initEul(3), 'XYZ')'; % todo: sim depends on q1 rn
q1 = q1(2:end);
x1 = [10; 0; 0; 0; 0; 0; q1; 0; 0; 0]; % todo: sim depends on x1 rn
u1_ = [0; 0; m*g; 0];

q2_ = angle2quat(0, 0, 0, 'XYZ')';
q2_ = q2_(2:end);
x2_ = [15; 15; -15; 0; 0; 0; q2_; 0; 0; 0];
u2_ = [0; 0; m*g; 0];

q3_ = angle2quat(0, 0, 0, 'XYZ')';
q3_ = q3_(2:end);
x3_ = [1; 0; 0; 0; 0; 0; q3_; 0; 0; 0];
u3_ = [0; 0; m*g; 0];

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

x_crit_ = [x1, x2_, x3_];
u_crit_ = [u1_, u2_, u3_];

discretize_interval = 0.0025; % loop time in seconds

%SIMULATE and create trajectory
fprintf("Creating Trajectory\n");
[x_set, u_set, t_set, Kset, tSegs, startTime, stopTime] = get_trajectory(10000, 15, x_crit_, u_crit_, [m; l; g], MOI, inputLimits, throttleConsts);

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
