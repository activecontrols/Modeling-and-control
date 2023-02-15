function [A_matrix,B_matrix] = get_A_B(x_current, u_current)
g = 9.81;
m = 5;
l = 0.5;
I1 = 10;
I5 = 100;
I9 = 100;
I2 = 0;
I3 = 0;
I4 = 0;
I6 = 0;
I7 = 0;
I8 = 0;

x_1 = x_current(1);
x_2 = x_current(2);
x_3 = x_current(3);
phi = x_current(4);
theta = x_current(5);
psi = x_current(6);
x_dot_1 = x_current(7);
x_dot_2 = x_current(8);
x_dot_3 = x_current(9);
phi_dot = x_current(10);
theta_dot = x_current(11);
psi_dot = x_current(12);

beta = u_current(1);
gamma = u_current(2);
F_T = u_current(3);
tau_RW = u_current(4);

A_matrix = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, g*cos(psi)*sin(theta), g*cos(theta)*sin(psi), 0, psi_dot, -theta_dot, 0, -x_dot_3, x_dot_2; 0, 0, 0, 0, -g*sin(psi)*sin(theta), g*cos(psi)*cos(theta), -psi_dot, 0, phi_dot, x_dot_3, 0, -x_dot_1; 0, 0, 0, 0, -g*cos(theta), 0, theta_dot, -phi_dot, 0, -x_dot_2, x_dot_1, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
B_matrix = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; -(F_T*cos(gamma)*sin(beta))/m, -(F_T*cos(beta)*sin(gamma))/m, (cos(beta)*cos(gamma))/m, 0; -(F_T*cos(beta))/m, 0, -sin(beta)/m, 0; 0, (F_T*cos(gamma))/m, sin(gamma)/m, 0; (F_T*l*cos(beta)*(I2*I6 - I3*I5))/(I1*I5*I9 - I1*I6*I8 - I2*I4*I9 + I2*I6*I7 + I3*I4*I8 - I3*I5*I7), -(F_T*l*cos(gamma)*(I2*I9 - I3*I8))/(I1*I5*I9 - I1*I6*I8 - I2*I4*I9 + I2*I6*I7 + I3*I4*I8 - I3*I5*I7), (l*sin(beta)*(I2*I6 - I3*I5))/(I1*I5*I9 - I1*I6*I8 - I2*I4*I9 + I2*I6*I7 + I3*I4*I8 - I3*I5*I7) - (l*sin(gamma)*(I2*I9 - I3*I8))/(I1*I5*I9 - I1*I6*I8 - I2*I4*I9 + I2*I6*I7 + I3*I4*I8 - I3*I5*I7), (I5*I9 - I6*I8)/(I1*I5*I9 - I1*I6*I8 - I2*I4*I9 + I2*I6*I7 + I3*I4*I8 - I3*I5*I7); -(F_T*l*cos(beta)*(I1*I6 - I3*I4))/(I1*I5*I9 - I1*I6*I8 - I2*I4*I9 + I2*I6*I7 + I3*I4*I8 - I3*I5*I7), (F_T*l*cos(gamma)*(I1*I9 - I3*I7))/(I1*I5*I9 - I1*I6*I8 - I2*I4*I9 + I2*I6*I7 + I3*I4*I8 - I3*I5*I7), (l*sin(gamma)*(I1*I9 - I3*I7))/(I1*I5*I9 - I1*I6*I8 - I2*I4*I9 + I2*I6*I7 + I3*I4*I8 - I3*I5*I7) - (l*sin(beta)*(I1*I6 - I3*I4))/(I1*I5*I9 - I1*I6*I8 - I2*I4*I9 + I2*I6*I7 + I3*I4*I8 - I3*I5*I7), -(I4*I9 - I6*I7)/(I1*I5*I9 - I1*I6*I8 - I2*I4*I9 + I2*I6*I7 + I3*I4*I8 - I3*I5*I7); (F_T*l*cos(beta)*(I1*I5 - I2*I4))/(I1*I5*I9 - I1*I6*I8 - I2*I4*I9 + I2*I6*I7 + I3*I4*I8 - I3*I5*I7), -(F_T*l*cos(gamma)*(I1*I8 - I2*I7))/(I1*I5*I9 - I1*I6*I8 - I2*I4*I9 + I2*I6*I7 + I3*I4*I8 - I3*I5*I7), (l*sin(beta)*(I1*I5 - I2*I4))/(I1*I5*I9 - I1*I6*I8 - I2*I4*I9 + I2*I6*I7 + I3*I4*I8 - I3*I5*I7) - (l*sin(gamma)*(I1*I8 - I2*I7))/(I1*I5*I9 - I1*I6*I8 - I2*I4*I9 + I2*I6*I7 + I3*I4*I8 - I3*I5*I7), (I4*I8 - I5*I7)/(I1*I5*I9 - I1*I6*I8 - I2*I4*I9 + I2*I6*I7 + I3*I4*I8 - I3*I5*I7)];

end

