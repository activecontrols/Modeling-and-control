function [A_matrix,B_matrix] = get_A_B(x_current)
syms x_1 x_2 x_3 phi theta psi x_dot_1 x_dot_2 x_dot_3 phi_dot theta_dot psi_dot beta gamma tau_RW
syms t g F_T l m
syms I [1 9]

x = [x_1; x_2; x_3; phi; theta; psi; x_dot_1; x_dot_2; x_dot_3; phi_dot; theta_dot; psi_dot];
u = [beta; gamma; F_T; tau_RW];
% deltaX = x - x_set;

x_dot_pos = [x_dot_1; x_dot_2; x_dot_3];
omega = [phi_dot; theta_dot; psi_dot];
x_ddot_pos = (1/m) * (-m*g*[cos(theta)*cos(psi); -cos(theta)*sin(psi); sin(theta)] + F_T*[cos(beta)*cos(gamma); -sin(beta); sin(gamma)])- cross(omega, x_dot_pos);
omega_dot = (inv([I1 I2 I3; I4 I5 I6; I7 I8 I9]))*(cross([-l; 0; 0], F_T*[cos(beta)*cos(gamma); -sin(beta); sin(gamma)]) + tau_RW*[1; 0; 0]);
x_dot = [x_dot_pos; omega; x_ddot_pos; omega_dot];

Jx = jacobian(x_dot, x);
Ju = jacobian(x_dot, u);

m = 5;
l = 0.5;
g = 9.81;
t = 0;
A = subs(Jx, x, xf);
A = subs(A);
B = subs(Ju, u, u0);
B = subs(B);
B = subs(B, [I1 I2 I3 I4 I5 I6 I7 I8 I9], [10 0 0 0 100 0 0 0 100]);

A_matrix = double(A);

B_matrix = double(vpa(B));
end

