function [x_set, u_set, t_set, K_set, tSegs, startTime, stopTime] = get_trajectory(numPoints, tf, xcritset, ucritset, constants, MOI, limits, throttleConsts)
    %Number of critical points in trajectory
    numCrits = size(xcritset, 2);
    %Full timeline initializations
    t_set = zeros(1, numCrits*numPoints);
    
    %initial time
    ti = 0;
    
    %time segment limits
    tSegs = linspace(ti, tf, numCrits);

    %maximum values for state and input used for Bryson's rule
    rmax = [.05; 1; 1];
    vmax = [.01; 1; 1];
    qmax = angle2quat(pi/120, pi/120, pi, 'ZYX')';
    qmax = qmax(2:end);
    omegamax = [pi/2; pi/2; pi/20];
    xmax = [rmax; vmax; qmax; omegamax];
    
    umax = [pi/12; pi/12; 100; 1];
    
    %Cost function Q and R matrices
    Q = diag(xmax.^-2);
    R = diag(umax.^-2);

    %trajectory set initialization
    xset = zeros(size(xcritset, 1), numPoints*numCrits);
    uset = zeros(size(ucritset, 1), numPoints*numCrits);
    K_set = zeros(size(uset, 1), size(xset, 1)+3, numCrits);

    %Create trajectories for different segments
    for i = 1:(numCrits - 1)
        if i > 1
            [xset(:, (1+numPoints*(i-1)):numPoints*i), uset(:, (1+numPoints*(i-1)):numPoints*i), t_set(1, (1+numPoints*(i-1)):numPoints*i), K_set(:, :, i)] = get_segment_traj(numPoints, tSegs(i), tSegs(i+1), xset(:, numPoints*(i-1)), xcritset(:, i+1), ucritset(:, i+1), Q, R, constants, MOI, limits, throttleConsts);
        else
            [xset(:, (1+numPoints*(i-1)):numPoints*i), uset(:, (1+numPoints*(i-1)):numPoints*i), t_set(1, (1+numPoints*(i-1)):numPoints*i), K_set(:, :, i)] = get_segment_traj(numPoints, tSegs(i), tSegs(i+1), xcritset(:, i), xcritset(:, i+1), ucritset(:, i+1), Q, R, constants, MOI, limits, throttleConsts);
        end
        xset(:, (1+numPoints*(i-1)):numPoints*i) = xset(:, (1+numPoints*(i-1)):numPoints*i) + xcritset(:, i+1);
    end
    
    %Get rid of duplicate points in trajectory
    [t_set, ia, ~] = unique(t_set, 'stable');
    x_set = xset(:, ia);
    u_set = uset(:, ia);
    startTime = t_set(1);
    stopTime = t_set(end);
end