function [xsegment, usegment, tsegment, Ksegment] = get_segment_traj(numPoints, ti, tf, xcrit1, xcrit2, ucrit2, Q, R, constants, MOI, limits, throttleConsts)

[x, u, xdot, Jx, Ju, consts, Jmat] = EOMS(throttleConsts);
J1 = Jmat(1, 1);
J2 = Jmat(1, 2);
J3 = Jmat(1, 3);
J4 = Jmat(2, 1);
J5 = Jmat(2, 2);
J6 = Jmat(2, 3);
J7 = Jmat(3, 1);
J8 = Jmat(3, 2);
J9 = Jmat(3, 3);

%Calculate A matrix in state space representation
A = subs(Jx, [x; u; consts], [zeros(length(xcrit2), 1); ucrit2; constants]);
A = double(subs(A, [J1 J2 J3; J4 J5 J6; J7 J8 J9], MOI));

%Calculate B matrix in state space representation
B = subs(Ju, [x; u; consts], [zeros(length(xcrit2), 1); ucrit2; constants]);
B = double(subs(B, [J1 J2 J3; J4 J5 J6; J7 J8 J9], MOI));

%Optimal control gain matrix K, solution S, and poles P
try
    [Ksegment, ~, ~] = lqr(A, B, Q, R);
catch e
    disp(e.message);
    error("LQR gain generation threw the error above!");
end

[xsegment, usegment, tsegment] = simulate(ti, tf, numPoints, Ksegment, constants, MOI, xcrit1-xcrit2, limits);

C = [eye(3), zeros(3, length(A)-3)];
Atilde = [A, zeros(size(A, 1), 3); C, zeros(size(C, 1), 3)];
Btilde = [B; zeros(3, size(B, 2))];
Qtilde = diag([diag(Q); (.01*ones(3, 1)).^-2]);
Rtilde = R; 

[Ksegment, ~, ~] = lqr(Atilde, Btilde, Qtilde, Rtilde);

end

function [xPlot, uPlot, tsegment] = simulate(ti, tf, numPoints, K, consts, MOI, xcrit, limits)
    %SIMULATE
    tsegment = linspace(ti, tf, numPoints);
    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    [~, xPlot] = ode45(@(t, x) deriv(t, x, K, consts, MOI, limits), tsegment, xcrit, opts);
    xPlot = xPlot';
    
    %Create simulated control inputs
    uPlot = zeros(size(K, 1), size(tsegment, 1));
    for i = 1:size(xPlot, 2)
        uPlot(:, i) = constrain(-K*xPlot(:, i), limits);
    end
end

function xdot = deriv(~, x, K, consts, MOI, limits)
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

    u = constrain(-K * x, limits);
    beta = u(1);
    gamma = u(2);
    throttle = u(3);
    tau_RW = u(4);

    xdot = [v1; v2; v3; -(g*m - sin(beta)*(2*q1*q3 + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q2)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) - cos(beta)*cos(gamma)*(2*(q2)^2 + 2*(q3)^2 - 1)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) + cos(beta)*sin(gamma)*(2*q1*q2 - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q3)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250))/m; (sin(beta)*(2*q2*q3 - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q1)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) + cos(beta)*sin(gamma)*(2*(q1)^2 + 2*(q3)^2 - 1)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) - cos(beta)*cos(gamma)*(2*q1*q2 + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q3)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250))/m; -(sin(beta)*(2*(q1)^2 + 2*(q2)^2 - 1)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) + cos(beta)*cos(gamma)*(2*q1*q3 - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q2)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) + cos(beta)*sin(gamma)*(2*q2*q3 + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q1)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250))/m; (omega1*sqrt(- q1^2 - q2^2 - q3^2 + 1))/2 - (omega2*q3)/2 + (omega3*q2)/2; (omega2*sqrt(- q1^2 - q2^2 - q3^2 + 1))/2 + (omega1*q3)/2 - (omega3*q1)/2; (omega3*sqrt(- q1^2 - q2^2 - q3^2 + 1))/2 - (omega1*q2)/2 + (omega2*q1)/2; (62500*J2^2*J6*omega2^2 + 62500*J2^2*J9*omega2*omega3 - 62500*J2*J3*J5*omega2^2 + 62500*J2*J3*J6*omega2*omega3 - 62500*J2*J3*J8*omega2*omega3 + 62500*J2*J3*J9*omega3^2 - 62500*J2*J5*J6*omega1*omega2 - 62500*J2*J6^2*omega1*omega3 - 62500*J4*J2*J6*omega1^2 + 62500*J1*J2*J6*omega1*omega2 + 7*l*cos(beta)*sin(gamma)*J2*J6*throttle^3 - 1295*l*cos(beta)*sin(gamma)*J2*J6*throttle^2 + 16750*l*cos(beta)*sin(gamma)*J2*J6*throttle - 62500*J2*J8*J9*omega1*omega2 - 62500*J2*J9^2*omega1*omega3 - 62500*J7*J2*J9*omega1^2 + 62500*J1*J2*J9*omega1*omega3 - 7*l*sin(beta)*J2*J9*throttle^3 + 1295*l*sin(beta)*J2*J9*throttle^2 - 16750*l*sin(beta)*J2*J9*throttle - 62500*J3^2*J5*omega2*omega3 - 62500*J3^2*J8*omega3^2 + 62500*J3*J5^2*omega1*omega2 + 62500*J3*J5*J6*omega1*omega3 + 62500*J4*J3*J5*omega1^2 - 62500*J1*J3*J5*omega1*omega2 - 7*l*cos(beta)*sin(gamma)*J3*J5*throttle^3 + 1295*l*cos(beta)*sin(gamma)*J3*J5*throttle^2 - 16750*l*cos(beta)*sin(gamma)*J3*J5*throttle + 62500*J3*J8^2*omega1*omega2 + 62500*J3*J8*J9*omega1*omega3 + 62500*J7*J3*J8*omega1^2 - 62500*J1*J3*J8*omega1*omega3 + 7*l*sin(beta)*J3*J8*throttle^3 - 1295*l*sin(beta)*J3*J8*throttle^2 + 16750*l*sin(beta)*J3*J8*throttle + 62500*J5^2*J9*omega2*omega3 - 62500*J5*J6*J8*omega2*omega3 + 62500*J5*J6*J9*omega3^2 - 62500*J5*J8*J9*omega2^2 - 62500*J5*J9^2*omega2*omega3 - 62500*J7*J5*J9*omega1*omega2 + 62500*J4*J5*J9*omega1*omega3 + 62500*tau_RW*J5*J9 - 62500*J6^2*J8*omega3^2 + 62500*J6*J8^2*omega2^2 + 62500*J6*J8*J9*omega2*omega3 + 62500*J7*J6*J8*omega1*omega2 - 62500*J4*J6*J8*omega1*omega3 - 62500*tau_RW*J6*J8)/(62500*(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7)); -(62500*J1^2*J6*omega1*omega2 + 62500*J1^2*J9*omega1*omega3 - 62500*J1*J3*J4*omega1*omega2 + 62500*J1*J3*J6*omega2*omega3 - 62500*J1*J3*J7*omega1*omega3 + 62500*J1*J3*J9*omega3^2 - 62500*J1*J4*J6*omega1^2 - 62500*J1*J6^2*omega1*omega3 - 62500*J5*J1*J6*omega1*omega2 + 62500*J2*J1*J6*omega2^2 + 7*l*cos(beta)*sin(gamma)*J1*J6*throttle^3 - 1295*l*cos(beta)*sin(gamma)*J1*J6*throttle^2 + 16750*l*cos(beta)*sin(gamma)*J1*J6*throttle - 62500*J1*J7*J9*omega1^2 - 62500*J1*J9^2*omega1*omega3 - 62500*J8*J1*J9*omega1*omega2 + 62500*J2*J1*J9*omega2*omega3 - 7*l*sin(beta)*J1*J9*throttle^3 + 1295*l*sin(beta)*J1*J9*throttle^2 - 16750*l*sin(beta)*J1*J9*throttle - 62500*J3^2*J4*omega2*omega3 - 62500*J3^2*J7*omega3^2 + 62500*J3*J4^2*omega1^2 + 62500*J3*J4*J6*omega1*omega3 + 62500*J5*J3*J4*omega1*omega2 - 62500*J2*J3*J4*omega2^2 - 7*l*cos(beta)*sin(gamma)*J3*J4*throttle^3 + 1295*l*cos(beta)*sin(gamma)*J3*J4*throttle^2 - 16750*l*cos(beta)*sin(gamma)*J3*J4*throttle + 62500*J3*J7^2*omega1^2 + 62500*J3*J7*J9*omega1*omega3 + 62500*J8*J3*J7*omega1*omega2 - 62500*J2*J3*J7*omega2*omega3 + 7*l*sin(beta)*J3*J7*throttle^3 - 1295*l*sin(beta)*J3*J7*throttle^2 + 16750*l*sin(beta)*J3*J7*throttle + 62500*J4^2*J9*omega1*omega3 - 62500*J4*J6*J7*omega1*omega3 + 62500*J4*J6*J9*omega3^2 - 62500*J4*J7*J9*omega1*omega2 - 62500*J4*J9^2*omega2*omega3 - 62500*J8*J4*J9*omega2^2 + 62500*J5*J4*J9*omega2*omega3 + 62500*tau_RW*J4*J9 - 62500*J6^2*J7*omega3^2 + 62500*J6*J7^2*omega1*omega2 + 62500*J6*J7*J9*omega2*omega3 + 62500*J8*J6*J7*omega2^2 - 62500*J5*J6*J7*omega2*omega3 - 62500*tau_RW*J6*J7)/(62500*(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7)); (62500*J1^2*J5*omega1*omega2 + 62500*J1^2*J8*omega1*omega3 - 62500*J1*J2*J4*omega1*omega2 + 62500*J1*J2*J5*omega2^2 - 62500*J1*J2*J7*omega1*omega3 + 62500*J1*J2*J8*omega2*omega3 - 62500*J1*J4*J5*omega1^2 - 62500*J1*J5^2*omega1*omega2 - 62500*J6*J1*J5*omega1*omega3 + 62500*J3*J1*J5*omega2*omega3 + 7*l*cos(beta)*sin(gamma)*J1*J5*throttle^3 - 1295*l*cos(beta)*sin(gamma)*J1*J5*throttle^2 + 16750*l*cos(beta)*sin(gamma)*J1*J5*throttle - 62500*J1*J7*J8*omega1^2 - 62500*J1*J8^2*omega1*omega2 - 62500*J9*J1*J8*omega1*omega3 + 62500*J3*J1*J8*omega3^2 - 7*l*sin(beta)*J1*J8*throttle^3 + 1295*l*sin(beta)*J1*J8*throttle^2 - 16750*l*sin(beta)*J1*J8*throttle - 62500*J2^2*J4*omega2^2 - 62500*J2^2*J7*omega2*omega3 + 62500*J2*J4^2*omega1^2 + 62500*J2*J4*J5*omega1*omega2 + 62500*J6*J2*J4*omega1*omega3 - 62500*J3*J2*J4*omega2*omega3 - 7*l*cos(beta)*sin(gamma)*J2*J4*throttle^3 + 1295*l*cos(beta)*sin(gamma)*J2*J4*throttle^2 - 16750*l*cos(beta)*sin(gamma)*J2*J4*throttle + 62500*J2*J7^2*omega1^2 + 62500*J2*J7*J8*omega1*omega2 + 62500*J9*J2*J7*omega1*omega3 - 62500*J3*J2*J7*omega3^2 + 7*l*sin(beta)*J2*J7*throttle^3 - 1295*l*sin(beta)*J2*J7*throttle^2 + 16750*l*sin(beta)*J2*J7*throttle + 62500*J4^2*J8*omega1*omega3 - 62500*J4*J5*J7*omega1*omega3 + 62500*J4*J5*J8*omega2*omega3 - 62500*J4*J7*J8*omega1*omega2 - 62500*J4*J8^2*omega2^2 - 62500*J9*J4*J8*omega2*omega3 + 62500*J6*J4*J8*omega3^2 + 62500*tau_RW*J4*J8 - 62500*J5^2*J7*omega2*omega3 + 62500*J5*J7^2*omega1*omega2 + 62500*J5*J7*J8*omega2^2 + 62500*J9*J5*J7*omega2*omega3 - 62500*J6*J5*J7*omega3^2 - 62500*tau_RW*J5*J7)/(62500*(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7))];
end

function unew = constrain(u, limits)
    unew = zeros(size(u, 1), size(u, 2));
    
    beta_min = limits(1, 1);
    beta_max = limits(1, 2);
    if u(1) > beta_max
        unew(1) = beta_max;
    elseif u(1) < beta_min
        unew(1) = beta_min;
    else
        unew(1) = u(1);
    end
    
    gamma_min = limits(2, 1);
    gamma_max = limits(2, 2);
    if u(2) > gamma_max
        unew(2) = gamma_max;
    elseif u(2) < gamma_min
        unew(2) = gamma_min;
    else
        unew(2) = u(2);
    end

    T_min = limits(3, 1);
    T_max = limits(3, 2);

    if u(3) > T_max
        unew(3) = T_max;
    elseif u(3) < T_min
        unew(3) = T_min;
    else
        unew(3) = u(3);
    end
    
    tau_min = limits(4, 1);
    tau_max = limits(4, 2);
    if u(4) > tau_max
        unew(4) = tau_max;
    elseif u(4) < tau_min
        unew(4) = tau_min;
    else
        unew(4) = u(4);
    end
end