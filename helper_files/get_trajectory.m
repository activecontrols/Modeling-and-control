function [xset, uset, tset, Kset, tSegs] = get_trajectory(numPoints, xcritset, ucritset, constants, MOI)
    mass = constants(1);
    grav = constants(3);
    numCrits = size(xcritset, 2);
    tset = zeros(1, numCrits*numPoints);
    
    ti = 0;
    tf = 5;
    tSegs = linspace(ti, tf, numCrits);

    rmax = [1; 1; 1];
    vmax = [1; 1; 1];
    qmax = angle2quat(pi/12, pi/12, pi/12, 'ZYX')';
    qmax = qmax(2:end);
    omegamax = [pi/2; pi/2; pi/2];
    xmax = [rmax; vmax; qmax; omegamax];
    
    umax = [pi/12; pi/12; mass*grav*5; 1];
    
    %Cost function Q and R matrices
    Q = diag(xmax.^-2);
    R = diag(umax.^-2);

    xset = zeros(size(xcritset, 1), numPoints*numCrits);
    uset = zeros(size(ucritset, 1), numPoints*numCrits);
    Kset = zeros(size(uset, 1), size(xset, 1)+3, numCrits);

    for i = 1:(numCrits - 1)
        [xset(:, (1+numPoints*(i-1)):numPoints*i), uset(:, (1+numPoints*(i-1)):numPoints*i), tset(1, (1+numPoints*(i-1)):numPoints*i), Kset(:, :, i)] = get_segment_traj(numPoints, tSegs(i), tSegs(i+1), xcritset(:, i), xcritset(:, i+1), ucritset(:, i+1), Q, R, constants, MOI);
        xset(:, (1+numPoints*(i-1)):numPoints*i) = xset(:, (1+numPoints*(i-1)):numPoints*i) + xcritset(:, i+1);
    end
end