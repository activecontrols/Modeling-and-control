clear;
clc;
numPoints = 10000;
filename = "Ellis_H50.csv";
raw = readmatrix(filename);
raw_t = raw(:, 1);
raw_thrust = raw(:, 2);
burn_time_vector = linspace(raw_t(1), raw_t(end), numPoints);
thrust = interp1(raw_t, raw_thrust, burn_time_vector);
syms t time a v r;
[raw_thrust_force, L] = thrust_series(raw_t, raw_thrust, 18);
Ftoft = subs(raw_thrust_force, time, t);

start_time = 0;
end_time = 100;
burn_start = end_time - L - 1;
time_vector = [linspace(start_time, burn_start, numPoints / 3), linspace(burn_start, burn_start + L, numPoints / 3), linspace(burn_start + L, end_time, numPoints / 3)];

m = 3;
v0 = 0;
r0 = 250;


a(t) = (Ftoft / m) - 9.81;
v(t) = int(a(t));
v = v - feval(v, 0) + v0;
r(t) = int(v(t));
r = r - feval(r, 0) + r0;

hold on;
subplot(1, 3, 1);
fplot(a(t), [time_vector(1), time_vector(end)]);
subplot(1, 3, 2);
fplot(v(t), [time_vector(1), time_vector(end)]);
subplot(1, 3, 3);
fplot(r(t), [time_vector(1), time_vector(end)]);