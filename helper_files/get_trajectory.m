function [x_set, u_set, t_set, K_set, tSegs, startTime, stopTime] = get_trajectory(paramArray)
    xcritset = paramArray{5};
    ucritset = paramArray{6};
    numPoints = paramArray{9};
    ti = paramArray{10};
    tf = paramArray{11};

    %Number of critical points in trajectory
    numCrits = size(xcritset, 2);

    %Full timeline initializations
    t_set = zeros(1, numCrits*numPoints);
    
    %time segment limits
    tSegs = linspace(ti, tf, numCrits);

    %trajectory set initialization
    xset = zeros(size(xcritset, 1), numPoints*numCrits);
    uset = zeros(size(ucritset, 1), numPoints*numCrits);
    K_set = zeros(size(uset, 1), size(xset, 1)+3, numCrits);

    %initialize segment specific parameter array
    segArray = paramArray;

    %Create trajectories for different segments
    for i = 1:(numCrits - 1)
        a = (1+numPoints*(i-1)); %bottom index
        b = numPoints*i; %top index

        %update segment speicific parameter array
        segArray{10} = tSegs(i);
        segArray{11} = tSegs(i+1);
        segArray{6} = ucritset(:, i+1);
        if i > 1
            segArray{5} = [xset(:, numPoints*(i-1)), xcritset(:, i+1)];
        else
            segArray{5} = [xcritset(:, i), xcritset(:, i+1)];
        end
        
        %call get_segment_traj.m
        [xset(:, a:b), uset(:, a:b), t_set(1, a:b), K_set(:, :, i)] = get_segment_traj(segArray);

        xset(:, a:b) = xset(:, a:b) + xcritset(:, i+1);
    end
    
    %Get rid of duplicate points in trajectory
    [t_set, ia, ~] = unique(t_set, 'stable');
    x_set = xset(:, ia);
    u_set = uset(:, ia);
    startTime = t_set(1);
    stopTime = t_set(end);
end