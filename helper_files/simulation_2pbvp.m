%Simulates and returns trajectory
function [xPlot, uPlot, tsegment, odeStopped] = simulation_2pbvp(ti, tf, numPoints, K, consts, MOI, xcrit, limits)
    %SIMULATE
    stopTime = 2.5;
    odeStopped = false;
    tsegment = linspace(ti, tf, numPoints);
    timer = tic;
    
    

    try

    catch

    end
end

function ydot = bvp_ode(~, y, K, consts, MOI, limits)
    % Extract constants
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
    
    % Extract individual states
    r1 = y(1);
    r2 = y(2);
    r3 = y(3);
    v1 = y(4);
    v2 = y(5);
    v3 = y(6);
    q1 = y(7);
    q2 = y(8);
    q3 = y(9);
    omega1 = x(10);
    omega2 = x(11);
    omega3 = x(12);

    % Extract costates
    p1 = y(13);
    p2 = y(14);
    p3 = y(15);
    p4 = y(16);
    p5 = y(17);
    p6 = y(18);
    p7 = y(19);
    p8 = y(20);
    p9 = y(21);
    p10 = y(22);
    p11 = y(23);
    p12 = y(24);

    % Extract control inputs
    % CALC U
    beta = u(1);
    gamma = u(2);
    throttle = u(3);
    tau_RW = u(4);



    xdot = [v1; 
            v2; 
            v3; 
            -(g*m - sin(beta)*(2*q1*q3 + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q2)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) - cos(beta)*cos(gamma)*(2*(q2)^2 + 2*(q3)^2 - 1)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) + cos(beta)*sin(gamma)*(2*q1*q2 - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q3)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250))/m; 
            (sin(beta)*(2*q2*q3 - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q1)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) + cos(beta)*sin(gamma)*(2*(q1)^2 + 2*(q3)^2 - 1)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) - cos(beta)*cos(gamma)*(2*q1*q2 + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q3)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250))/m; 
            -(sin(beta)*(2*(q1)^2 + 2*(q2)^2 - 1)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) + cos(beta)*cos(gamma)*(2*q1*q3 - 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q2)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250) + cos(beta)*sin(gamma)*(2*q2*q3 + 2*conj(sqrt(- q1^2 - q2^2 - q3^2 + 1))*q1)*((7*throttle^3)/62500 - (259*throttle^2)/12500 + (67*throttle)/250))/m; 
            (omega1*sqrt(- q1^2 - q2^2 - q3^2 + 1))/2 - (omega2*q3)/2 + (omega3*q2)/2; 
            (omega2*sqrt(- q1^2 - q2^2 - q3^2 + 1))/2 + (omega1*q3)/2 - (omega3*q1)/2; 
            (omega3*sqrt(- q1^2 - q2^2 - q3^2 + 1))/2 - (omega1*q2)/2 + (omega2*q1)/2; 
            (62500*J2^2*J6*omega2^2 + 62500*J2^2*J9*omega2*omega3 - 62500*J2*J3*J5*omega2^2 + 62500*J2*J3*J6*omega2*omega3 - 62500*J2*J3*J8*omega2*omega3 + 62500*J2*J3*J9*omega3^2 - 62500*J2*J5*J6*omega1*omega2 - 62500*J2*J6^2*omega1*omega3 - 62500*J4*J2*J6*omega1^2 + 62500*J1*J2*J6*omega1*omega2 + 7*l*cos(beta)*sin(gamma)*J2*J6*throttle^3 - 1295*l*cos(beta)*sin(gamma)*J2*J6*throttle^2 + 16750*l*cos(beta)*sin(gamma)*J2*J6*throttle - 62500*J2*J8*J9*omega1*omega2 - 62500*J2*J9^2*omega1*omega3 - 62500*J7*J2*J9*omega1^2 + 62500*J1*J2*J9*omega1*omega3 - 7*l*sin(beta)*J2*J9*throttle^3 + 1295*l*sin(beta)*J2*J9*throttle^2 - 16750*l*sin(beta)*J2*J9*throttle - 62500*J3^2*J5*omega2*omega3 - 62500*J3^2*J8*omega3^2 + 62500*J3*J5^2*omega1*omega2 + 62500*J3*J5*J6*omega1*omega3 + 62500*J4*J3*J5*omega1^2 - 62500*J1*J3*J5*omega1*omega2 - 7*l*cos(beta)*sin(gamma)*J3*J5*throttle^3 + 1295*l*cos(beta)*sin(gamma)*J3*J5*throttle^2 - 16750*l*cos(beta)*sin(gamma)*J3*J5*throttle + 62500*J3*J8^2*omega1*omega2 + 62500*J3*J8*J9*omega1*omega3 + 62500*J7*J3*J8*omega1^2 - 62500*J1*J3*J8*omega1*omega3 + 7*l*sin(beta)*J3*J8*throttle^3 - 1295*l*sin(beta)*J3*J8*throttle^2 + 16750*l*sin(beta)*J3*J8*throttle + 62500*J5^2*J9*omega2*omega3 - 62500*J5*J6*J8*omega2*omega3 + 62500*J5*J6*J9*omega3^2 - 62500*J5*J8*J9*omega2^2 - 62500*J5*J9^2*omega2*omega3 - 62500*J7*J5*J9*omega1*omega2 + 62500*J4*J5*J9*omega1*omega3 + 62500*tau_RW*J5*J9 - 62500*J6^2*J8*omega3^2 + 62500*J6*J8^2*omega2^2 + 62500*J6*J8*J9*omega2*omega3 + 62500*J7*J6*J8*omega1*omega2 - 62500*J4*J6*J8*omega1*omega3 - 62500*tau_RW*J6*J8)/(62500*(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7)); 
            -(62500*J1^2*J6*omega1*omega2 + 62500*J1^2*J9*omega1*omega3 - 62500*J1*J3*J4*omega1*omega2 + 62500*J1*J3*J6*omega2*omega3 - 62500*J1*J3*J7*omega1*omega3 + 62500*J1*J3*J9*omega3^2 - 62500*J1*J4*J6*omega1^2 - 62500*J1*J6^2*omega1*omega3 - 62500*J5*J1*J6*omega1*omega2 + 62500*J2*J1*J6*omega2^2 + 7*l*cos(beta)*sin(gamma)*J1*J6*throttle^3 - 1295*l*cos(beta)*sin(gamma)*J1*J6*throttle^2 + 16750*l*cos(beta)*sin(gamma)*J1*J6*throttle - 62500*J1*J7*J9*omega1^2 - 62500*J1*J9^2*omega1*omega3 - 62500*J8*J1*J9*omega1*omega2 + 62500*J2*J1*J9*omega2*omega3 - 7*l*sin(beta)*J1*J9*throttle^3 + 1295*l*sin(beta)*J1*J9*throttle^2 - 16750*l*sin(beta)*J1*J9*throttle - 62500*J3^2*J4*omega2*omega3 - 62500*J3^2*J7*omega3^2 + 62500*J3*J4^2*omega1^2 + 62500*J3*J4*J6*omega1*omega3 + 62500*J5*J3*J4*omega1*omega2 - 62500*J2*J3*J4*omega2^2 - 7*l*cos(beta)*sin(gamma)*J3*J4*throttle^3 + 1295*l*cos(beta)*sin(gamma)*J3*J4*throttle^2 - 16750*l*cos(beta)*sin(gamma)*J3*J4*throttle + 62500*J3*J7^2*omega1^2 + 62500*J3*J7*J9*omega1*omega3 + 62500*J8*J3*J7*omega1*omega2 - 62500*J2*J3*J7*omega2*omega3 + 7*l*sin(beta)*J3*J7*throttle^3 - 1295*l*sin(beta)*J3*J7*throttle^2 + 16750*l*sin(beta)*J3*J7*throttle + 62500*J4^2*J9*omega1*omega3 - 62500*J4*J6*J7*omega1*omega3 + 62500*J4*J6*J9*omega3^2 - 62500*J4*J7*J9*omega1*omega2 - 62500*J4*J9^2*omega2*omega3 - 62500*J8*J4*J9*omega2^2 + 62500*J5*J4*J9*omega2*omega3 + 62500*tau_RW*J4*J9 - 62500*J6^2*J7*omega3^2 + 62500*J6*J7^2*omega1*omega2 + 62500*J6*J7*J9*omega2*omega3 + 62500*J8*J6*J7*omega2^2 - 62500*J5*J6*J7*omega2*omega3 - 62500*tau_RW*J6*J7)/(62500*(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7)); 
            (62500*J1^2*J5*omega1*omega2 + 62500*J1^2*J8*omega1*omega3 - 62500*J1*J2*J4*omega1*omega2 + 62500*J1*J2*J5*omega2^2 - 62500*J1*J2*J7*omega1*omega3 + 62500*J1*J2*J8*omega2*omega3 - 62500*J1*J4*J5*omega1^2 - 62500*J1*J5^2*omega1*omega2 - 62500*J6*J1*J5*omega1*omega3 + 62500*J3*J1*J5*omega2*omega3 + 7*l*cos(beta)*sin(gamma)*J1*J5*throttle^3 - 1295*l*cos(beta)*sin(gamma)*J1*J5*throttle^2 + 16750*l*cos(beta)*sin(gamma)*J1*J5*throttle - 62500*J1*J7*J8*omega1^2 - 62500*J1*J8^2*omega1*omega2 - 62500*J9*J1*J8*omega1*omega3 + 62500*J3*J1*J8*omega3^2 - 7*l*sin(beta)*J1*J8*throttle^3 + 1295*l*sin(beta)*J1*J8*throttle^2 - 16750*l*sin(beta)*J1*J8*throttle - 62500*J2^2*J4*omega2^2 - 62500*J2^2*J7*omega2*omega3 + 62500*J2*J4^2*omega1^2 + 62500*J2*J4*J5*omega1*omega2 + 62500*J6*J2*J4*omega1*omega3 - 62500*J3*J2*J4*omega2*omega3 - 7*l*cos(beta)*sin(gamma)*J2*J4*throttle^3 + 1295*l*cos(beta)*sin(gamma)*J2*J4*throttle^2 - 16750*l*cos(beta)*sin(gamma)*J2*J4*throttle + 62500*J2*J7^2*omega1^2 + 62500*J2*J7*J8*omega1*omega2 + 62500*J9*J2*J7*omega1*omega3 - 62500*J3*J2*J7*omega3^2 + 7*l*sin(beta)*J2*J7*throttle^3 - 1295*l*sin(beta)*J2*J7*throttle^2 + 16750*l*sin(beta)*J2*J7*throttle + 62500*J4^2*J8*omega1*omega3 - 62500*J4*J5*J7*omega1*omega3 + 62500*J4*J5*J8*omega2*omega3 - 62500*J4*J7*J8*omega1*omega2 - 62500*J4*J8^2*omega2^2 - 62500*J9*J4*J8*omega2*omega3 + 62500*J6*J4*J8*omega3^2 + 62500*tau_RW*J4*J8 - 62500*J5^2*J7*omega2*omega3 + 62500*J5*J7^2*omega1*omega2 + 62500*J5*J7*J8*omega2^2 + 62500*J9*J5*J7*omega2*omega3 - 62500*J6*J5*J7*omega3^2 - 62500*tau_RW*J5*J7)/(62500*(J1*J5*J9 - J1*J6*J8 - J2*J4*J9 + J2*J6*J7 + J3*J4*J8 - J3*J5*J7))];
end

function res = bcfun(x, xcrit)
    res = xcrit - x;
end

function [value,isterminal,direction] = time_EventsFcn(~, ~, timer, stopTime) 
    value = 1; % The value that we want to be zero 
    if stopTime - toc(timer) < 0 % Halt if he has not finished in 3     
        error("ODE45:runtimeEvent", "Integration stopped: time longer than %f seconds", stopTime)
    end 
    isterminal = 1;  % Halt integration  
    direction = 0; % The zero can be approached from either direction 
end