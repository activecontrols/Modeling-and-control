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

%initial and final time
ti = 0;
tf = 27;

%number of points in trajectory
numPoints = 10000;

%genetic algorithm settings
genOn = false;
popSize = 250;
mut_rate1 = 0.20;
mut_rate2 = 0.01;
gen_cut = 0.8;
elite_cut = 0;
costFunc_weights = [0.1667, 0.1667]; % [r overshoot, r final error]

%Enumerate operating points for trajectory below:
initEul = [pi/2, pi/6, pi/6];
q1 = angle2quat(initEul(1), initEul(2), initEul(3), 'XYZ')';
q1 = q1(2:end);
x1 = [10; 10; 10; 0; 0; 0; q1; 0; 0; 0];
u1 = [0; 0; m*g; 0];

q2 = angle2quat(0, 0, 0, 'XYZ')';
q2 = q2(2:end);
x2 = [5; 5; 0; 0; 0; 0; q2; 0; 0; 0];
u2 = [0; 0; m*g; 0];

q3 = angle2quat(0, 0, 0, 'XYZ')';
q3 = q3(2:end);
x3 = [1; 0; 0; 0; 0; 0; q3; 0; 0; 0];
u3 = [0; 0; m*g; 0];

x_crits = [x1, x2, x3];
u_crits = [u1, u2, u3];

%Maximum and minimum values for input for simulation
inputLimits = [-pi/12, pi/12; 
               -pi/12, pi/12; 
                    0, 100;
                   -2, 2];

%[min, max] values for state for simulation
stateLimits = [-100, 100;
                  0, 100; 
               -100, 100; 
                NaN, NaN;
                NaN, NaN;
                NaN, NaN; 
                NaN, NaN; 
                NaN, NaN; 
                NaN, NaN; 
                NaN, NaN; 
                NaN, NaN; 
                NaN, NaN];

%Throttle -> Thrust force equation constants
throttleConsts = [-0.000112; 0.02072; -.268];

%Delay for servo angle inputs
betaInputDelay = 0.001;
gammaInputDelay = 0.001;

%Determine diagonal entries of cost function Q and R matrices using brysons rule (seed solution set of GA)
rmax = [1; 1; 1];
vmax = [.1; 1; 1];
qmax = angle2quat(pi/1200, pi/1200, pi, 'ZYX')';
qmax = qmax(2:end);
omegamax = [pi/20; pi/20; pi/20];
xmax = [rmax; vmax; qmax; omegamax];

umax = [pi/12; pi/12; 1000; 1];

Qbry = (xmax.^-2);
Rbry = (umax.^-2);

% Create cell array to store params
paramArray = {MOI, l, m, g, x_crits, u_crits, inputLimits, stateLimits, throttleConsts, numPoints, ti, tf, Qbry, Rbry, genOn, popSize, mut_rate1, mut_rate2, gen_cut, elite_cut, costFunc_weights};

%SIMULATE and create trajectory
begTrajGen = tic;
fprintf("Creating Trajectory\n");
[x_set, u_set, t_set, Kset, tSegs, startTime, stopTime] = get_trajectory(paramArray);
fprintf("Time to get trajectory = %f s", toc(begTrajGen))

%plot trajectory
plotTrajectory(x_set, u_set, t_set, 0.75, 50, [m; l; g]);

%plot states
figure(2)
plot(t_set, x_set(1,:), '-', 'color', "#0072BD")
hold on
plot(t_set, x_set(2,:), '-', 'color', "#D95319")
hold on
plot(t_set, x_set(3,:), '-', 'color', "#EDB120")
hold on
plot(t_set, x_set(1,:), '--', 'color', "#0072BD")
hold on
plot(t_set, x_set(2,:), '--', 'color', "#D95319")
hold on
plot(t_set, x_set(3,:), '--', 'color', "#EDB120")
hold on
plot(t_set, x_set(7,:), '-', 'color', "#7E2F8E")
hold on
plot(t_set, x_set(8,:), '-', 'color', "#77AC30")
hold on
plot(t_set, x_set(9,:), '-', 'color', "#A2142F")
hold on
plot(t_set, x_set(7,:), '--', 'color', "#7E2F8E")
hold on
plot(t_set, x_set(8,:), '--', 'color', "#77AC30")
hold on
plot(t_set, x_set(9,:), '--', 'color', "#A2142F")
hold on
title('States vs Time')
xlabel('Time (s)')
ylabel('States')
legend('r1', 'r2', 'r3', 'v1', 'v2', 'v3', 'q1', 'q2', 'q3', 'w1', 'w2', 'w3')
grid on

fprintf("Done initializing!\n");
cd ../