grav = 9.81;
acc = 30;
burnTime = 5;

syms time
linAcc = [piecewise(time<=burnTime, acc - grav, time>burnTime, -grav); 0; 0];
angAcc = [0; 0; 0];
timeVector = linspace(0, 20, 200);
xInit = [0.5; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

traj = calc_trajectory(timeVector, xInit, linAcc, angAcc);

syms x_1 x_2 x_3 phi theta psi x_dot_1 x_dot_2 x_dot_3 phi_dot theta_dot psi_dot beta gamma tau_RW
syms t g F_T l m
syms I [1 9]

x = [x_1; x_2; x_3; phi; theta; psi; x_dot_1; x_dot_2; x_dot_3; phi_dot; theta_dot; psi_dot];
u = [beta; gamma; F_T; tau_RW];


x_dot_pos = [x_dot_1; x_dot_2; x_dot_3];
omega = [phi_dot; theta_dot; psi_dot];
x_ddot_pos = (1/m) * (-m*g*[cos(theta)*cos(psi); -cos(theta)*sin(psi); sin(theta)] + F_T*[cos(beta)*cos(gamma); -sin(beta); sin(gamma)])- cross(omega, x_dot_pos);
omega_dot = (inv([I1 I2 I3; I4 I5 I6; I7 I8 I9]))*(cross([-l; 0; 0], F_T*[cos(beta)*cos(gamma); -sin(beta); sin(gamma)]) + tau_RW*[1; 0; 0]);
x_dot = [x_dot_pos; omega; x_ddot_pos; omega_dot];

Jx = jacobian(x_dot, x);
Ju = jacobian(x_dot, u);
g = 9.81;
m = 3;
F_T = m*linAcc(1);
l = 0.5;
sumList = zeros(1, length(traj));
for i = 1:length(traj)
    x_current = traj(i, :)';
    A = subs(Jx, x, x_current);
    A = subs(A);
    sum = 0;
    for j = 1:length(traj)
        x_iter = traj(j, :)';
        J_current = subs(Jx, x, x_iter);
        J_current = subs(J_current);
        sum = sum + statSize(A - J_current);
    end
    sumList(i) = sum;
    fprintf("A#: %d\tSum: %d\n", i, sum);
end

plot(sumList);

function sum = statSize(M)
sum = 0;
    for i = size(M, 1)
        for j = size(M, 2)
            sum = sum + M(i, j)^2;
        end
    end
sum;
end