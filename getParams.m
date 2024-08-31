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
% % initEul = [pi/3, pi/12, pi/12];
% initEul = [0, 0, 0];
% % order is roll, pitch, yaw?
% q1 = angle2quat(initEul(1), initEul(2), initEul(3), 'XYZ')'; % todo: sim depends on q1 rn
% q1 = q1(2:end);
% x1 = [10; 0; 0; 0; 0; 0; q1; 0; 0; 0]; % todo: sim depends on x1 rn
% u1 = [0; 0; m*g; 0];
% 
% q2 = angle2quat(0, 0, 0, 'XYZ')';
% q2 = q2(2:end);
% x2 = [15; 15; -15; 0; 0; 0; q2; 0; 0; 0];
% u2 = [0; 0; m*g; 0];
% 
% q3 = angle2quat(0, 0, 0, 'XYZ')';
% q3 = q3(2:end);
% x3 = [1; 0; 0; 0; 0; 0; q3; 0; 0; 0];
% u3 = [0; 0; m*g; 0];

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
% inputLimits = [-pi/12, pi/12; -pi/12, pi/12; 0, 100; -2, 2];
% no reaction wheel input:
inputLimits = [-pi/12, pi/12; -pi/12, pi/12; 0, 100; -0.1, 0.1];

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

discretize_interval = .0025; % loop time in seconds

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



%%
% 
% Z = [0.623613796917501, 0.6247204231021626, 0.6250942374047459, 0.6269209757431385, 0.6031725419329573, 0.6026457057512153, 0.6022818154413538, 0.5989020879350337, 0.6323166673564201, 0.6311795970749097, 0.6323718405835453, 0.6318482433605124, 0.6030302597928479, 0.6023771411257819, 0.6001707239445706, 0.6021167452393531]
% 
% plot(Z, ones(size(Z)), 'r.')