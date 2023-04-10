clear;

cd ./PDM_integration/
raw = readmatrix("/parameters.xlsx");
cd ../
MOI = raw(1:3, 2:4);
l = raw(4, 2);
m = raw(5, 2);

g = 9.81;

syms x1 x2 x3 v1 v2 v3 q1 q2 q3 omega1 omega2 omega3 betaAng gammaAng F_T tau_RW
syms J [1 9]
r = [x1; x2; x3];
v = [v1; v2; v3];
q0 = sqrt(1 - q1^2 - q2^2 - q3^2);
qVec = [q0; q1; q2; q3];
qChange = [q1; q2; q3];
omegaB = [omega1; omega2; omega3];
x = [r; v; qChange; omegaB];
u = [betaAng; gammaAng; F_T; tau_RW];

Jmat = [J1 J2 J3; J4 J5 J6; J7 J8 J9];

TB = F_T * [cos(betaAng)*cos(gammaAng); cos(betaAng)*sin(gammaAng); -sin(betaAng)];
uSym = [betaAng; gammaAng; F_T; tau_RW];
xSym = [r; v; qChange; omegaB];
CIB = quatRot(qVec);
CBI = CIB';

rdot = v;
FI = CBI*TB;
vdot = (1 / m) * (FI + [-m*g; 0; 0]);
qdot = 0.5 * bigOmega(omegaB)*qVec;
qdot = qdot(2:end);
rTB = [-l; 0; 0];
MB = zetaCross(rTB)*TB + [tau_RW; 0; 0];
omegaBdot = Jmat \ (MB - zetaCross(omegaB)*Jmat*omegaB);

fSym = [rdot; vdot; qdot; omegaBdot];

max_x = [1; 1; 1; 1; 1; 1; .01; .01; .01; .1; .1; .1];
% max_x = ones(1, 12);
max_u = [pi/12; pi/12; .01; .001];
% max_u = ones(1, 4);
Q = zeros(length(x), length(x));
R = zeros(length(u), length(u));
for i = 1:length(x)
    Q(i, i) = 1/(max_x(i)^2);
end
for j = 1:length(u)
    R(j, j) = 1/(max_u(j)^2); 
end

qi = angle2quat(0, 0, 0, 'ZYX')';
qi = qi(2:end);

qf = angle2quat(0, 0, 0, 'ZYX')';
qf = qf(2:end);

x0 = [100; 0; 0; 0; 0; 0; qi; 0; 0; 0];
xf = [l; 0; 0; 0; 0; 0; qf; 0; 0; 0];
u0 = [0; 0; 0; 0];
uf = [0; 0; m*g; 0];
endTime = 10;
%% 
fprintf("Creating Trajectory\n");
[xref, uref, tRef] = get_trajectory(x0, xf, u0, uf, endTime, 1000, m, l, MOI, g);
% plot(uref(1, :));

fprintf("Done initializing!\n");

%% 
[A_stateError, B_stateError] = AB(xf, uf, m, MOI, l, g);
Atilde = [A_stateError, zeros(length(A_stateError)); -1*eye(length(A_stateError)), zeros(length(A_stateError))];
Btilde = [B_stateError; zeros(length(A_stateError), size(B_stateError, 2))];
Qtilde = diag([max_x.^(-2); zeros(length(max_x), 1)]);
Rtilde = diag([max_u.^(-2); zeros(length(max_u)-length(max_x), 1)]);
ucontrol = length(Atilde) - rank(ctrb(Atilde, Btilde))
[Ktilde, ~, ~] = lqr(Atilde, Btilde, Qtilde, Rtilde);

%% 
function CIB = quatRot(q)
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    CIB = [1-2*(q2^2 + q3^2), 2*(q1*q2 + q0*q3), 2*(q1*q3 - q0*q2);
        2*(q1*q2 - q0*q3), 1-2*(q1^2 + q3^2), 2*(q2*q3 + q0*q1);
        2*(q1*q3 + q0*q2), 2*(q2*q3 - q0*q1), 1-2*(q1^2 + q2^2)];
end

function zc = zetaCross(zeta)
    zc = [0 -zeta(3) zeta(2); zeta(3) 0 -zeta(1); -zeta(2) zeta(1) 0];
end
function bigO = bigOmega(zeta)
    bigO = [0 -zeta(1) -zeta(2) -zeta(3); zeta(1) 0 zeta(3) -zeta(2); zeta(2) -zeta(3) 0 zeta(1); zeta(3) zeta(2) -zeta(1) 0];
end
function xdot = f(t, x, u, m, J, l, g)
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    v1 = x(4);
    v2 = x(5);
    v3 = x(6);
    q1 = x(7);
    q2 = x(8);
    q3 = x(9);
    omega1 = x(10);
    omega2 = x(11);
    omega3 = x(12);

    betaAng = u(1);
    gammaAng = u(2);
    F_T = u(3);
    tau_RW = u(4);

    J1 = J(1);
    J2 = J(2);
    J3 = J(3);
    J4 = J(4);
    J5 = J(5);
    J6 = J(6);
    J7 = J(7);
    J8 = J(8);
    J9 = J(9);

    xdot = [v1; v2; v3; -(g*m + F_T*sin(betaAng)*(2*conj(q1)*conj(q3) + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q2)) - F_T*cos(betaAng)*sin(gammaAng)*(2*conj(q1)*conj(q2) - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q3)) + F_T*cos(betaAng)*cos(gammaAng)*(2*(conj(q2))^2 + 2*(conj(q3))^2 - 1))/m; -(F_T*sin(betaAng)*(2*conj(q2)*conj(q3) - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q1)) - F_T*cos(betaAng)*cos(gammaAng)*(2*conj(q1)*conj(q2) + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q3)) + F_T*cos(betaAng)*sin(gammaAng)*(2*(conj(q1))^2 + 2*(conj(q3))^2 - 1))/m; (F_T*sin(betaAng)*(2*(conj(q1))^2 + 2*(conj(q2))^2 - 1) + F_T*cos(betaAng)*cos(gammaAng)*(2*conj(q1)*conj(q3) - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q2)) + F_T*cos(betaAng)*sin(gammaAng)*(2*conj(q2)*conj(q3) + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q1)))/m; (omega1*sqrt(- q1^2 - q2^2 - q3^2 + 1))/2 - (omega2*q3)/2 + (omega3*q2)/2; (omega2*sqrt(- q1^2 - q2^2 - q3^2 + 1))/2 + (omega1*q3)/2 - (omega3*q1)/2; (omega3*sqrt(- q1^2 - q2^2 - q3^2 + 1))/2 - (omega1*q2)/2 + (omega2*q1)/2; (J2^2*J6*omega2^2 + J2^2*J9*omega2*omega3 - J2*J3*J5*omega2^2 + J2*J3*J6*omega2*omega3 - J2*J3*J8*omega2*omega3 + J2*J3*J9*omega3^2 - J2*J5*J6*omega1*omega2 - J2*J6^2*omega1*omega3 - J4*J2*J6*omega1^2 + J1*J2*J6*omega1*omega2 - F_T*l*cos(betaAng)*sin(gammaAng)*J2*J6 - J2*J8*J9*omega1*omega2 - J2*J9^2*omega1*omega3 - J7*J2*J9*omega1^2 + J1*J2*J9*omega1*omega3 + F_T*l*sin(betaAng)*J2*J9 - J3^2*J5*omega2*omega3 - J3^2*J8*omega3^2 + J3*J5^2*omega1*omega2 + J3*J5*J6*omega1*omega3 + J4*J3*J5*omega1^2 - J1*J3*J5*omega1*omega2 + F_T*l*cos(betaAng)*sin(gammaAng)*J3*J5 + J3*J8^2*omega1*omega2 + J3*J8*J9*omega1*omega3 + J7*J3*J8*omega1^2 - J1*J3*J8*omega1*omega3 - F_T*l*sin(betaAng)*J3*J8 + J5^2*J9*omega2*omega3 - J5*J6*J8*omega2*omega3 + J5*J6*J9*omega3^2 - J5*J8*J9*omega2^2 - J5*J9^2*omega2*omega3 - J7*J5*J9*omega1*omega2 + J4*J5*J9*omega1*omega3 + tau_RW*J5*J9 - J6^2*J8*omega3^2 + J6*J8^2*omega2^2 + J6*J8*J9*omega2*omega3 + J7*J6*J8*omega1*omega2 - J4*J6*J8*omega1*omega3 - tau_RW*J6*J8)/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7); -(J1^2*J6*omega1*omega2 + J1^2*J9*omega1*omega3 - J1*J3*J4*omega1*omega2 + J1*J3*J6*omega2*omega3 - J1*J3*J7*omega1*omega3 + J1*J3*J9*omega3^2 - J1*J4*J6*omega1^2 - J1*J6^2*omega1*omega3 - J5*J1*J6*omega1*omega2 + J2*J1*J6*omega2^2 - F_T*l*cos(betaAng)*sin(gammaAng)*J1*J6 - J1*J7*J9*omega1^2 - J1*J9^2*omega1*omega3 - J8*J1*J9*omega1*omega2 + J2*J1*J9*omega2*omega3 + F_T*l*sin(betaAng)*J1*J9 - J3^2*J4*omega2*omega3 - J3^2*J7*omega3^2 + J3*J4^2*omega1^2 + J3*J4*J6*omega1*omega3 + J5*J3*J4*omega1*omega2 - J2*J3*J4*omega2^2 + F_T*l*cos(betaAng)*sin(gammaAng)*J3*J4 + J3*J7^2*omega1^2 + J3*J7*J9*omega1*omega3 + J8*J3*J7*omega1*omega2 - J2*J3*J7*omega2*omega3 - F_T*l*sin(betaAng)*J3*J7 + J4^2*J9*omega1*omega3 - J4*J6*J7*omega1*omega3 + J4*J6*J9*omega3^2 - J4*J7*J9*omega1*omega2 - J4*J9^2*omega2*omega3 - J8*J4*J9*omega2^2 + J5*J4*J9*omega2*omega3 + tau_RW*J4*J9 - J6^2*J7*omega3^2 + J6*J7^2*omega1*omega2 + J6*J7*J9*omega2*omega3 + J8*J6*J7*omega2^2 - J5*J6*J7*omega2*omega3 - tau_RW*J6*J7)/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7); (J1^2*J5*omega1*omega2 + J1^2*J8*omega1*omega3 - J1*J2*J4*omega1*omega2 + J1*J2*J5*omega2^2 - J1*J2*J7*omega1*omega3 + J1*J2*J8*omega2*omega3 - J1*J4*J5*omega1^2 - J1*J5^2*omega1*omega2 - J6*J1*J5*omega1*omega3 + J3*J1*J5*omega2*omega3 - F_T*l*cos(betaAng)*sin(gammaAng)*J1*J5 - J1*J7*J8*omega1^2 - J1*J8^2*omega1*omega2 - J9*J1*J8*omega1*omega3 + J3*J1*J8*omega3^2 + F_T*l*sin(betaAng)*J1*J8 - J2^2*J4*omega2^2 - J2^2*J7*omega2*omega3 + J2*J4^2*omega1^2 + J2*J4*J5*omega1*omega2 + J6*J2*J4*omega1*omega3 - J3*J2*J4*omega2*omega3 + F_T*l*cos(betaAng)*sin(gammaAng)*J2*J4 + J2*J7^2*omega1^2 + J2*J7*J8*omega1*omega2 + J9*J2*J7*omega1*omega3 - J3*J2*J7*omega3^2 - F_T*l*sin(betaAng)*J2*J7 + J4^2*J8*omega1*omega3 - J4*J5*J7*omega1*omega3 + J4*J5*J8*omega2*omega3 - J4*J7*J8*omega1*omega2 - J4*J8^2*omega2^2 - J9*J4*J8*omega2*omega3 + J6*J4*J8*omega3^2 + tau_RW*J4*J8 - J5^2*J7*omega2*omega3 + J5*J7^2*omega1*omega2 + J5*J7*J8*omega2^2 + J9*J5*J7*omega2*omega3 - J6*J5*J7*omega3^2 - tau_RW*J5*J7)/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7)];
end
function [A, B] = AB(x, u, m, J, l, g)
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    v1 = x(4);
    v2 = x(5);
    v3 = x(6);
    q1 = x(7);
    q2 = x(8);
    q3 = x(9);
    omega1 = x(10);
    omega2 = x(11);
    omega3 = x(12);

    betaAng = u(1);
    gammaAng = u(2);
    F_T = u(3);

    J1 = J(1, 1);
    J2 = J(1, 2);
    J3 = J(1, 3);
    J4 = J(2, 1);
    J5 = J(2, 2);
    J6 = J(2, 3);
    J7 = J(3, 1);
    J8 = J(3, 2);
    J9 = J(3, 3);

    A = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -(F_T*sin(betaAng)*(2*conj(q3) - (2*conj(q1)*conj(q2))/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))) - F_T*cos(betaAng)*sin(gammaAng)*(2*conj(q2) + (2*conj(q1)*conj(q3))/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))))/m, -(F_T*sin(betaAng)*(2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1)) - (2*(conj(q2))^2)/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))) + 4*F_T*cos(betaAng)*cos(gammaAng)*conj(q2) - F_T*cos(betaAng)*sin(gammaAng)*(2*conj(q1) + (2*conj(q2)*conj(q3))/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))))/m, -(F_T*sin(betaAng)*(2*conj(q1) - (2*conj(q2)*conj(q3))/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))) + 4*F_T*cos(betaAng)*cos(gammaAng)*conj(q3) + F_T*cos(betaAng)*sin(gammaAng)*(2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1)) - (2*(conj(q3))^2)/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))))/m, 0, 0, 0; 0, 0, 0, 0, 0, 0, (F_T*sin(betaAng)*(2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1)) - (2*(conj(q1))^2)/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))) - 4*F_T*cos(betaAng)*conj(q1)*sin(gammaAng) + F_T*cos(betaAng)*cos(gammaAng)*(2*conj(q2) - (2*conj(q1)*conj(q3))/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))))/m, -(F_T*sin(betaAng)*(2*conj(q3) + (2*conj(q1)*conj(q2))/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))) - F_T*cos(betaAng)*cos(gammaAng)*(2*conj(q1) - (2*conj(q2)*conj(q3))/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))))/m, -(F_T*sin(betaAng)*(2*conj(q2) + (2*conj(q1)*conj(q3))/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))) + 4*F_T*cos(betaAng)*conj(q3)*sin(gammaAng) - F_T*cos(betaAng)*cos(gammaAng)*(2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1)) - (2*(conj(q3))^2)/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))))/m, 0, 0, 0; 0, 0, 0, 0, 0, 0, (4*F_T*conj(q1)*sin(betaAng) + F_T*cos(betaAng)*sin(gammaAng)*(2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1)) - (2*(conj(q1))^2)/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))) + F_T*cos(betaAng)*cos(gammaAng)*(2*conj(q3) + (2*conj(q1)*conj(q2))/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))))/m, (4*F_T*conj(q2)*sin(betaAng) - F_T*cos(betaAng)*cos(gammaAng)*(2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1)) - (2*(conj(q2))^2)/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))) + F_T*cos(betaAng)*sin(gammaAng)*(2*conj(q3) - (2*conj(q1)*conj(q2))/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))))/m, (F_T*cos(betaAng)*cos(gammaAng)*(2*conj(q1) + (2*conj(q2)*conj(q3))/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))) + F_T*cos(betaAng)*sin(gammaAng)*(2*conj(q2) - (2*conj(q1)*conj(q3))/conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))))/m, 0, 0, 0; 0, 0, 0, 0, 0, 0, -(omega1*q1)/(2*sqrt(- q1^2 - q2^2 - q3^2 + 1)), omega3/2 - (omega1*q2)/(2*sqrt(- q1^2 - q2^2 - q3^2 + 1)), - omega2/2 - (omega1*q3)/(2*sqrt(- q1^2 - q2^2 - q3^2 + 1)), sqrt(- q1^2 - q2^2 - q3^2 + 1)/2, -q3/2, q2/2; 0, 0, 0, 0, 0, 0, - omega3/2 - (omega2*q1)/(2*sqrt(- q1^2 - q2^2 - q3^2 + 1)), -(omega2*q2)/(2*sqrt(- q1^2 - q2^2 - q3^2 + 1)), omega1/2 - (omega2*q3)/(2*sqrt(- q1^2 - q2^2 - q3^2 + 1)), q3/2, sqrt(- q1^2 - q2^2 - q3^2 + 1)/2, -q1/2; 0, 0, 0, 0, 0, 0, omega2/2 - (omega3*q1)/(2*sqrt(- q1^2 - q2^2 - q3^2 + 1)), - omega1/2 - (omega3*q2)/(2*sqrt(- q1^2 - q2^2 - q3^2 + 1)), -(omega3*q3)/(2*sqrt(- q1^2 - q2^2 - q3^2 + 1)), -q2/2, q1/2, sqrt(- q1^2 - q2^2 - q3^2 + 1)/2; 0, 0, 0, 0, 0, 0, 0, 0, 0, (J3*J5^2*omega2 - J2*J6^2*omega3 + J3*J8^2*omega2 - J2*J9^2*omega3 + J1*J2*J6*omega2 - J1*J3*J5*omega2 - 2*J2*J4*J6*omega1 + 2*J3*J4*J5*omega1 + J1*J2*J9*omega3 - J1*J3*J8*omega3 - J2*J5*J6*omega2 + J3*J5*J6*omega3 - 2*J2*J7*J9*omega1 + 2*J3*J7*J8*omega1 - J2*J8*J9*omega2 + J4*J5*J9*omega3 - J4*J6*J8*omega3 + J3*J8*J9*omega3 - J5*J7*J9*omega2 + J6*J7*J8*omega2)/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7), (omega3*J2^2*J9 + 2*J6*omega2*J2^2 - 2*omega2*J2*J3*J5 - omega3*J2*J3*J8 + J6*omega3*J2*J3 - J6*omega1*J2*J5 - omega1*J2*J8*J9 + J1*J6*omega1*J2 - omega3*J3^2*J5 + omega1*J3*J5^2 - J1*omega1*J3*J5 + omega1*J3*J8^2 + omega3*J5^2*J9 - 2*omega2*J5*J8*J9 - J6*omega3*J5*J8 - omega3*J5*J9^2 - J7*omega1*J5*J9 + 2*J6*omega2*J8^2 + J6*omega3*J8*J9 + J6*J7*omega1*J8)/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7), -(- omega2*J2^2*J9 - omega2*J2*J3*J6 - 2*omega3*J2*J3*J9 + J8*omega2*J2*J3 + omega1*J2*J6^2 + omega1*J2*J9^2 - J1*omega1*J2*J9 + omega2*J3^2*J5 + 2*J8*omega3*J3^2 - omega1*J3*J5*J6 - J8*omega1*J3*J9 + J1*J8*omega1*J3 - omega2*J5^2*J9 - 2*omega3*J5*J6*J9 + J8*omega2*J5*J6 + omega2*J5*J9^2 - J4*omega1*J5*J9 + 2*J8*omega3*J6^2 - J8*omega2*J6*J9 + J4*J8*omega1*J6)/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7); 0, 0, 0, 0, 0, 0, 0, 0, 0, -(omega2*J1^2*J6 + omega3*J1^2*J9 - 2*omega1*J1*J4*J6 - J3*omega2*J1*J4 - omega3*J1*J6^2 - J5*omega2*J1*J6 - 2*omega1*J1*J7*J9 - J3*omega3*J1*J7 - omega3*J1*J9^2 - J8*omega2*J1*J9 + omega3*J4^2*J9 + 2*J3*omega1*J4^2 - omega3*J4*J6*J7 + J3*omega3*J4*J6 - omega2*J4*J7*J9 + J3*J5*omega2*J4 + omega2*J6*J7^2 + 2*J3*omega1*J7^2 + J3*omega3*J7*J9 + J3*J8*omega2*J7)/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7), -(J1^2*J6*omega1 - J3^2*J4*omega3 + J6*J7^2*omega1 - J4*J9^2*omega3 - J1*J3*J4*omega1 + 2*J1*J2*J6*omega2 - 2*J2*J3*J4*omega2 + J1*J3*J6*omega3 - J1*J5*J6*omega1 + J3*J4*J5*omega1 + J1*J2*J9*omega3 - J2*J3*J7*omega3 - J1*J8*J9*omega1 + J3*J7*J8*omega1 + J4*J5*J9*omega3 - J4*J7*J9*omega1 - J5*J6*J7*omega3 - 2*J4*J8*J9*omega2 + 2*J6*J7*J8*omega2 + J6*J7*J9*omega3)/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7), (- omega1*J1^2*J9 - omega2*J1*J3*J6 - 2*omega3*J1*J3*J9 + J7*omega1*J1*J3 + omega1*J1*J6^2 + omega1*J1*J9^2 - J2*omega2*J1*J9 + omega2*J3^2*J4 + 2*J7*omega3*J3^2 - omega1*J3*J4*J6 - J7*omega1*J3*J9 + J2*J7*omega2*J3 - omega1*J4^2*J9 - 2*omega3*J4*J6*J9 + J7*omega1*J4*J6 + omega2*J4*J9^2 - J5*omega2*J4*J9 + 2*J7*omega3*J6^2 - J7*omega2*J6*J9 + J5*J7*omega2*J6)/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7); 0, 0, 0, 0, 0, 0, 0, 0, 0, (omega2*J1^2*J5 + omega3*J1^2*J8 - 2*omega1*J1*J4*J5 - J2*omega2*J1*J4 - omega2*J1*J5^2 - J6*omega3*J1*J5 - 2*omega1*J1*J7*J8 - J2*omega3*J1*J7 - omega2*J1*J8^2 - J9*omega3*J1*J8 + omega3*J4^2*J8 + 2*J2*omega1*J4^2 - omega3*J4*J5*J7 + J2*omega2*J4*J5 - omega2*J4*J7*J8 + J2*J6*omega3*J4 + omega2*J5*J7^2 + 2*J2*omega1*J7^2 + J2*omega2*J7*J8 + J2*J9*omega3*J7)/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7), -(- omega1*J1^2*J5 - 2*omega2*J1*J2*J5 - omega3*J1*J2*J8 + J4*omega1*J1*J2 + omega1*J1*J5^2 - J3*omega3*J1*J5 + omega1*J1*J8^2 + omega3*J2^2*J7 + 2*J4*omega2*J2^2 - J4*omega1*J2*J5 - omega1*J2*J7*J8 + J3*J4*omega3*J2 + omega3*J5^2*J7 - omega1*J5*J7^2 - 2*omega2*J5*J7*J8 - J9*omega3*J5*J7 - J4*omega3*J5*J8 + J4*omega1*J7*J8 + 2*J4*omega2*J8^2 + J4*J9*omega3*J8)/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7), (J1^2*J8*omega1 - J2^2*J7*omega2 + J4^2*J8*omega1 - J5^2*J7*omega2 - J1*J2*J7*omega1 + J1*J3*J5*omega2 - J2*J3*J4*omega2 + J1*J2*J8*omega2 - J1*J5*J6*omega1 + J2*J4*J6*omega1 + 2*J1*J3*J8*omega3 - 2*J2*J3*J7*omega3 - J4*J5*J7*omega1 - J1*J8*J9*omega1 + J2*J7*J9*omega1 + J4*J5*J8*omega2 + 2*J4*J6*J8*omega3 - 2*J5*J6*J7*omega3 - J4*J8*J9*omega2 + J5*J7*J9*omega2)/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7)];
    B = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; -(F_T*cos(betaAng)*(2*conj(q1)*conj(q3) + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q2)) + F_T*sin(betaAng)*sin(gammaAng)*(2*conj(q1)*conj(q2) - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q3)) - F_T*cos(gammaAng)*sin(betaAng)*(2*(conj(q2))^2 + 2*(conj(q3))^2 - 1))/m, (F_T*cos(betaAng)*cos(gammaAng)*(2*conj(q1)*conj(q2) - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q3)) + F_T*cos(betaAng)*sin(gammaAng)*(2*(conj(q2))^2 + 2*(conj(q3))^2 - 1))/m, -(sin(betaAng)*(2*conj(q1)*conj(q3) + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q2)) - cos(betaAng)*sin(gammaAng)*(2*conj(q1)*conj(q2) - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q3)) + cos(betaAng)*cos(gammaAng)*(2*(conj(q2))^2 + 2*(conj(q3))^2 - 1))/m, 0; -(F_T*cos(betaAng)*(2*conj(q2)*conj(q3) - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q1)) + F_T*cos(gammaAng)*sin(betaAng)*(2*conj(q1)*conj(q2) + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q3)) - F_T*sin(betaAng)*sin(gammaAng)*(2*(conj(q1))^2 + 2*(conj(q3))^2 - 1))/m, -(F_T*cos(betaAng)*sin(gammaAng)*(2*conj(q1)*conj(q2) + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q3)) + F_T*cos(betaAng)*cos(gammaAng)*(2*(conj(q1))^2 + 2*(conj(q3))^2 - 1))/m, -(sin(betaAng)*(2*conj(q2)*conj(q3) - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q1)) - cos(betaAng)*cos(gammaAng)*(2*conj(q1)*conj(q2) + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q3)) + cos(betaAng)*sin(gammaAng)*(2*(conj(q1))^2 + 2*(conj(q3))^2 - 1))/m, 0; -(F_T*cos(gammaAng)*sin(betaAng)*(2*conj(q1)*conj(q3) - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q2)) - F_T*cos(betaAng)*(2*(conj(q1))^2 + 2*(conj(q2))^2 - 1) + F_T*sin(betaAng)*sin(gammaAng)*(2*conj(q2)*conj(q3) + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q1)))/m, (F_T*cos(betaAng)*cos(gammaAng)*(2*conj(q2)*conj(q3) + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q1)) - F_T*cos(betaAng)*sin(gammaAng)*(2*conj(q1)*conj(q3) - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q2)))/m, (sin(betaAng)*(2*(conj(q1))^2 + 2*(conj(q2))^2 - 1) + cos(betaAng)*cos(gammaAng)*(2*conj(q1)*conj(q3) - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q2)) + cos(betaAng)*sin(gammaAng)*(2*conj(q2)*conj(q3) + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*conj(q1)))/m, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; (F_T*J2*J9*l*cos(betaAng) - F_T*J3*J8*l*cos(betaAng) + F_T*J2*J6*l*sin(betaAng)*sin(gammaAng) - F_T*J3*J5*l*sin(betaAng)*sin(gammaAng))/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7), -(F_T*J2*J6*l*cos(betaAng)*cos(gammaAng) - F_T*J3*J5*l*cos(betaAng)*cos(gammaAng))/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7), (J2*J9*l*sin(betaAng) - J3*J8*l*sin(betaAng) - J2*J6*l*cos(betaAng)*sin(gammaAng) + J3*J5*l*cos(betaAng)*sin(gammaAng))/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7), (J5*J9 - J6*J8)/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7); -(F_T*J1*J9*l*cos(betaAng) - F_T*J3*J7*l*cos(betaAng) + F_T*J1*J6*l*sin(betaAng)*sin(gammaAng) - F_T*J3*J4*l*sin(betaAng)*sin(gammaAng))/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7), (F_T*J1*J6*l*cos(betaAng)*cos(gammaAng) - F_T*J3*J4*l*cos(betaAng)*cos(gammaAng))/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7), -(J1*J9*l*sin(betaAng) - J3*J7*l*sin(betaAng) - J1*J6*l*cos(betaAng)*sin(gammaAng) + J3*J4*l*cos(betaAng)*sin(gammaAng))/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7), -(J4*J9 - J6*J7)/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7); (F_T*J1*J8*l*cos(betaAng) - F_T*J2*J7*l*cos(betaAng) + F_T*J1*J5*l*sin(betaAng)*sin(gammaAng) - F_T*J2*J4*l*sin(betaAng)*sin(gammaAng))/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7), -(F_T*J1*J5*l*cos(betaAng)*cos(gammaAng) - F_T*J2*J4*l*cos(betaAng)*cos(gammaAng))/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7), (J1*J8*l*sin(betaAng) - J2*J7*l*sin(betaAng) - J1*J5*l*cos(betaAng)*sin(gammaAng) + J2*J4*l*cos(betaAng)*sin(gammaAng))/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7), (J4*J8 - J5*J7)/(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7)];
end