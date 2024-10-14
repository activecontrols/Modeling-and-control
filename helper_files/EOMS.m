function [x, u, xdot, Jx, Ju, consts, Jmat, aug_eoms, u_sol] = EOMS(throttleConsts)
    %Declare symbolic variables
    syms r1 r2 r3 v1 v2 v3 q1 q2 q3 omega1 omega2 omega3 betaAng gammaAng throttle tau_RW m l g
    syms J [1 9]
    syms c [1 12] % costates
    syms R [1 16]% arbitrary weights for lagrangian L(x,u)

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
    p1 = throttleConsts(1);
    p2 = throttleConsts(2);
    p3 = throttleConsts(3);
    F_T = p1*throttle^3 + p2*throttle^2 + p3*throttle;
    
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

    %Augmented eoms: ydot = [xdot, cdot]' where cdot is the time derivative of the costate vector c
    % Pre-allocate variable sizes
    cdot = zeros(size(xdot));
    du = zeros(size(u));
    u_sol = du;
    
    c = [c1; c2; c3; c4; c5; c6; c7; c8; c9; c10; c11; c12];
    Rmat = [R1 R2 R3 R4; R5 R6 R7 R8; R9 R10 R11 R12; R13 R14 R15 R16];

    L = 0.5 * u' * Rmat * u; %Lagrangian
    H = L + c' * xdot; %Hamiltonian
    
    for i = 1:12
        cdot(i) = -diff(H, x(i));
    end
    for i = 1:4
        du(i) = diff(H, u(i));
        u_sol{i} = solve(du(i), u(i));
    end
    
    aug_eoms = [xdot; cdot];
    
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