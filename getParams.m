cd ./PDM_integration/
raw = readmatrix("/parameters.xlsx");
cd ../
MOI = raw(1:3, 2:4);
l = raw(4, 2);
m = raw(5, 2);

syms x_1 x_2 x_3 phi_angle_angle theta_angle psi_angle x_dot_1 x_dot_2 x_dot_3 phi_angle_dot theta_angle_dot psi_angle_dot tau_RW beta_angle gamma_angle;

syms t g F_T;
syms I [1 9];

x = [x_1; x_2; x_3; phi_angle_angle; theta_angle; psi_angle; x_dot_1; x_dot_2; x_dot_3; phi_angle_dot; theta_angle_dot; psi_angle_dot];
u = [beta_angle; gamma_angle; F_T; tau_RW];

x_dot_pos = [x_dot_1; x_dot_2; x_dot_3];
omega = [phi_angle_dot; theta_angle_dot; psi_angle_dot];
x_ddot_pos = (1/m) * (-m*g*[cos(theta_angle)*cos(psi_angle); -cos(theta_angle)*sin(psi_angle); sin(theta_angle)] + F_T*[cos(beta_angle)*cos(gamma_angle); -sin(beta_angle); sin(gamma_angle)])- cross(omega, x_dot_pos);
omega_dot = (inv([I1 I2 I3; I4 I5 I6; I7 I8 I9]))*(cross([-l; 0; 0], F_T*[cos(beta_angle)*cos(gamma_angle); -sin(beta_angle); sin(gamma_angle)]) + tau_RW*[1; 0; 0]);
x_dot = [x_dot_pos; omega; x_ddot_pos; omega_dot];

Jx = jacobian(x_dot, x);
Jx = vpa(Jx);
Ju = jacobian(x_dot, u);
Ju = vpa(Ju);

max_x = [.0001; .1; .1; pi/1200; pi/12; pi/12; .001; .1; .1; pi/12; pi/240; pi/240];
max_u = [pi/12; pi/12; .1; .5];
Q = zeros(length(x), length(x));
R = zeros(length(u), length(u));
for i = 1:length(x)
    Q(i, i) = 1/(max_x(i)^2);
end
for j = 1:length(u)
    R(j, j) = 1/(max_u(j)^2);
end

x0 = [60; 0; 0; 0; pi/120; 0; 0; 0; 0; 0; 0; 0];
xf = [l; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
endTime = 3;

g = 9.81;