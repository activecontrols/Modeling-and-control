function [x, u, xdot, Jx, Ju, consts, Jmat] = EOMS(throttleConsts)
    %Declare symbolic variables
    syms r1 r2 r3 v1 v2 v3 q1 q2 q3 omega1 omega2 omega3 betaAng gammaAng throttle tau_RW m l g
    syms J [1 9]
    
    %Create constants array
    consts = [m; l; g];
    
    %Earth-frame radius
    r = [r1; r2; r3];
    %Earth-frame velocity
    v = [v1; v2; v3];
    
    %Quaternion scalar relation formula
    q0 = sqrt(1 - q1^2 - q2^2 - q3^2);
    %Full quaternion vector
    qVec = [q0; q1; q2; q3];
    %Chopped Quaternion vector
    qChange = [q1; q2; q3];
    %Angular velocity (taken wrt E-frame, represented in B-frame)
    omegaB = [omega1; omega2; omega3];

    %Thrust force in terms of throttle input
    % p1 = throttleConsts(1);
    % p2 = throttleConsts(2);
    % p3 = throttleConsts(3);
    % F_T = p1*throttle^3 + p2*throttle^2 + p3*throttle;
    F_T = throttleConsts(1) ./ (1 + exp(-throttleConsts(2)*(throttle-throttleConsts(3))));
    
    %State vector
    x = [r; v; qChange; omegaB];
    %Input vector
    u = [betaAng; gammaAng; throttle; tau_RW];
    
    %Moment of inertia tensor
    Jmat = [J1 J2 J3; J4 J5 J6; J7 J8 J9];
    
    %Body-frame thrust vector
    TB = F_T * [cos(betaAng)*cos(gammaAng); cos(betaAng)*sin(gammaAng); -sin(betaAng)];
    
    %Direction-cosine rotation matrix (Earth to body)
    CIB = quatRot(qVec);
    %Direction-cosine rotation matrix (Body to earth)
    CBI = CIB';
    
    %Derivative of position = v
    rdot = v;
    
    %Net force in earth frame
    FI = CBI*TB - [m*g; 0; 0];
    
    %Derivative of velocity = acceleration = F / m
    vdot = FI / m;
    
    %Quaternion matrix
    bigO = [q0 -q3 q2; q3 q0 -q1; -q2 q1 q0];
    %Derivative of chopped quaternion array
    qdot = simplify(0.5 * bigO * omegaB);
    
    %Moment arm between center of mass and thrust point
    rTB = [-l; 0; 0];
    
    %Net moment in body frame
    MB = zetaCross(rTB)*TB + [tau_RW; 0; 0];
    
    %Derivative of angular velocity
    omegaBdot = Jmat \ (MB - zetaCross(omegaB)*Jmat*omegaB);
    
    %Symbolic state derivative
    xdot = [rdot; vdot; qdot; omegaBdot];
    
    %Jacobian of state derivative wrt state (Symbolic A matrix)
    Jx = jacobian(xdot, x);
    %Jacobian of state derivative wrt inputs (Symbolic B matrix)
    Ju = jacobian(xdot, u);
end

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

function xdot = deriv(~, x, K, consts, MOI)
    m = consts(1);
    l = consts(2);
    g = consts(3);
    J1 = MOI(1, 1);
    J2 = MOI(1, 2);
    J3 = MOI(1, 3);
    J4 = MOI(2, 1);
    J5 = MOI(2, 2);
    J6 = MOI(2, 3);
    J7 = MOI(3, 1);
    J8 = MOI(3, 2);
    J9 = MOI(3, 3);
    v1 = x(4);
    v2 = x(5);
    v3 = x(6);
    q1 = x(7);
    q2 = x(8);
    q3 = x(9);
    omega1 = x(10);
    omega2 = x(11);
    omega3 = x(12);

    u = constrain(-K * x);
    beta = u(1);
    gamma = u(2);
    throttle = u(3);
    tau_RW = u(4);

    xdot = [v1; v2; v3; -(g*m - sin(beta)*(2*q1*q3 + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q2)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) - cos(beta)*cos(gamma)*(2*(q2)^2 + 2*(q3)^2 - 1)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) + cos(beta)*sin(gamma)*(2*q1*q2 - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q3)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250))/m; (sin(beta)*(2*q2*q3 - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q1)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) + cos(beta)*sin(gamma)*(2*(q1)^2 + 2*(q3)^2 - 1)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) - cos(beta)*cos(gamma)*(2*q1*q2 + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q3)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250))/m; -(sin(beta)*(2*(q1)^2 + 2*(q2)^2 - 1)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) + cos(beta)*cos(gamma)*(2*q1*q3 - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q2)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) + cos(beta)*sin(gamma)*(2*q2*q3 + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q1)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250))/m; (omega1*sqrt(- q1^2 - q2^2 - q3^2 + 1))/2 - (omega2*q3)/2 + (omega3*q2)/2; (omega2*sqrt(- q1^2 - q2^2 - q3^2 + 1))/2 + (omega1*q3)/2 - (omega3*q1)/2; (omega3*sqrt(- q1^2 - q2^2 - q3^2 + 1))/2 - (omega1*q2)/2 + (omega2*q1)/2; (62500*J2^2*J6*omega2^2 + 62500*J2^2*J9*omega2*omega3 - 62500*J2*J3*J5*omega2^2 + 62500*J2*J3*J6*omega2*omega3 - 62500*J2*J3*J8*omega2*omega3 + 62500*J2*J3*J9*omega3^2 - 62500*J2*J5*J6*omega1*omega2 - 62500*J2*J6^2*omega1*omega3 - 62500*J4*J2*J6*omega1^2 + 62500*J1*J2*J6*omega1*omega2 + 7*l*cos(beta)*sin(gamma)*J2*J6*throttle^3 - 1295*l*cos(beta)*sin(gamma)*J2*J6*throttle^2 + 16750*l*cos(beta)*sin(gamma)*J2*J6*throttle - 62500*J2*J8*J9*omega1*omega2 - 62500*J2*J9^2*omega1*omega3 - 62500*J7*J2*J9*omega1^2 + 62500*J1*J2*J9*omega1*omega3 - 7*l*sin(beta)*J2*J9*throttle^3 + 1295*l*sin(beta)*J2*J9*throttle^2 - 16750*l*sin(beta)*J2*J9*throttle - 62500*J3^2*J5*omega2*omega3 - 62500*J3^2*J8*omega3^2 + 62500*J3*J5^2*omega1*omega2 + 62500*J3*J5*J6*omega1*omega3 + 62500*J4*J3*J5*omega1^2 - 62500*J1*J3*J5*omega1*omega2 - 7*l*cos(beta)*sin(gamma)*J3*J5*throttle^3 + 1295*l*cos(beta)*sin(gamma)*J3*J5*throttle^2 - 16750*l*cos(beta)*sin(gamma)*J3*J5*throttle + 62500*J3*J8^2*omega1*omega2 + 62500*J3*J8*J9*omega1*omega3 + 62500*J7*J3*J8*omega1^2 - 62500*J1*J3*J8*omega1*omega3 + 7*l*sin(beta)*J3*J8*throttle^3 - 1295*l*sin(beta)*J3*J8*throttle^2 + 16750*l*sin(beta)*J3*J8*throttle + 62500*J5^2*J9*omega2*omega3 - 62500*J5*J6*J8*omega2*omega3 + 62500*J5*J6*J9*omega3^2 - 62500*J5*J8*J9*omega2^2 - 62500*J5*J9^2*omega2*omega3 - 62500*J7*J5*J9*omega1*omega2 + 62500*J4*J5*J9*omega1*omega3 + 62500*tau_RW*J5*J9 - 62500*J6^2*J8*omega3^2 + 62500*J6*J8^2*omega2^2 + 62500*J6*J8*J9*omega2*omega3 + 62500*J7*J6*J8*omega1*omega2 - 62500*J4*J6*J8*omega1*omega3 - 62500*tau_RW*J6*J8)/(62500*(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7)); -(62500*J1^2*J6*omega1*omega2 + 62500*J1^2*J9*omega1*omega3 - 62500*J1*J3*J4*omega1*omega2 + 62500*J1*J3*J6*omega2*omega3 - 62500*J1*J3*J7*omega1*omega3 + 62500*J1*J3*J9*omega3^2 - 62500*J1*J4*J6*omega1^2 - 62500*J1*J6^2*omega1*omega3 - 62500*J5*J1*J6*omega1*omega2 + 62500*J2*J1*J6*omega2^2 + 7*l*cos(beta)*sin(gamma)*J1*J6*throttle^3 - 1295*l*cos(beta)*sin(gamma)*J1*J6*throttle^2 + 16750*l*cos(beta)*sin(gamma)*J1*J6*throttle - 62500*J1*J7*J9*omega1^2 - 62500*J1*J9^2*omega1*omega3 - 62500*J8*J1*J9*omega1*omega2 + 62500*J2*J1*J9*omega2*omega3 - 7*l*sin(beta)*J1*J9*throttle^3 + 1295*l*sin(beta)*J1*J9*throttle^2 - 16750*l*sin(beta)*J1*J9*throttle - 62500*J3^2*J4*omega2*omega3 - 62500*J3^2*J7*omega3^2 + 62500*J3*J4^2*omega1^2 + 62500*J3*J4*J6*omega1*omega3 + 62500*J5*J3*J4*omega1*omega2 - 62500*J2*J3*J4*omega2^2 - 7*l*cos(beta)*sin(gamma)*J3*J4*throttle^3 + 1295*l*cos(beta)*sin(gamma)*J3*J4*throttle^2 - 16750*l*cos(beta)*sin(gamma)*J3*J4*throttle + 62500*J3*J7^2*omega1^2 + 62500*J3*J7*J9*omega1*omega3 + 62500*J8*J3*J7*omega1*omega2 - 62500*J2*J3*J7*omega2*omega3 + 7*l*sin(beta)*J3*J7*throttle^3 - 1295*l*sin(beta)*J3*J7*throttle^2 + 16750*l*sin(beta)*J3*J7*throttle + 62500*J4^2*J9*omega1*omega3 - 62500*J4*J6*J7*omega1*omega3 + 62500*J4*J6*J9*omega3^2 - 62500*J4*J7*J9*omega1*omega2 - 62500*J4*J9^2*omega2*omega3 - 62500*J8*J4*J9*omega2^2 + 62500*J5*J4*J9*omega2*omega3 + 62500*tau_RW*J4*J9 - 62500*J6^2*J7*omega3^2 + 62500*J6*J7^2*omega1*omega2 + 62500*J6*J7*J9*omega2*omega3 + 62500*J8*J6*J7*omega2^2 - 62500*J5*J6*J7*omega2*omega3 - 62500*tau_RW*J6*J7)/(62500*(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7)); (62500*J1^2*J5*omega1*omega2 + 62500*J1^2*J8*omega1*omega3 - 62500*J1*J2*J4*omega1*omega2 + 62500*J1*J2*J5*omega2^2 - 62500*J1*J2*J7*omega1*omega3 + 62500*J1*J2*J8*omega2*omega3 - 62500*J1*J4*J5*omega1^2 - 62500*J1*J5^2*omega1*omega2 - 62500*J6*J1*J5*omega1*omega3 + 62500*J3*J1*J5*omega2*omega3 + 7*l*cos(beta)*sin(gamma)*J1*J5*throttle^3 - 1295*l*cos(beta)*sin(gamma)*J1*J5*throttle^2 + 16750*l*cos(beta)*sin(gamma)*J1*J5*throttle - 62500*J1*J7*J8*omega1^2 - 62500*J1*J8^2*omega1*omega2 - 62500*J9*J1*J8*omega1*omega3 + 62500*J3*J1*J8*omega3^2 - 7*l*sin(beta)*J1*J8*throttle^3 + 1295*l*sin(beta)*J1*J8*throttle^2 - 16750*l*sin(beta)*J1*J8*throttle - 62500*J2^2*J4*omega2^2 - 62500*J2^2*J7*omega2*omega3 + 62500*J2*J4^2*omega1^2 + 62500*J2*J4*J5*omega1*omega2 + 62500*J6*J2*J4*omega1*omega3 - 62500*J3*J2*J4*omega2*omega3 - 7*l*cos(beta)*sin(gamma)*J2*J4*throttle^3 + 1295*l*cos(beta)*sin(gamma)*J2*J4*throttle^2 - 16750*l*cos(beta)*sin(gamma)*J2*J4*throttle + 62500*J2*J7^2*omega1^2 + 62500*J2*J7*J8*omega1*omega2 + 62500*J9*J2*J7*omega1*omega3 - 62500*J3*J2*J7*omega3^2 + 7*l*sin(beta)*J2*J7*throttle^3 - 1295*l*sin(beta)*J2*J7*throttle^2 + 16750*l*sin(beta)*J2*J7*throttle + 62500*J4^2*J8*omega1*omega3 - 62500*J4*J5*J7*omega1*omega3 + 62500*J4*J5*J8*omega2*omega3 - 62500*J4*J7*J8*omega1*omega2 - 62500*J4*J8^2*omega2^2 - 62500*J9*J4*J8*omega2*omega3 + 62500*J6*J4*J8*omega3^2 + 62500*tau_RW*J4*J8 - 62500*J5^2*J7*omega2*omega3 + 62500*J5*J7^2*omega1*omega2 + 62500*J5*J7*J8*omega2^2 + 62500*J9*J5*J7*omega2*omega3 - 62500*J6*J5*J7*omega3^2 - 62500*tau_RW*J5*J7)/(62500*(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7))];
end