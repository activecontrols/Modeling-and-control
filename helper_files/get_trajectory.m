function x_set = get_trajectory(t)
    h0 = 10;
    g = 9.81;
    t1 = 1;
    b = 1962/1019;
    a = (1019*exp(1962/1019))/200;
    if (t < t1)
        x_set = [h0 - (g/2)*t^2; zeros(5, 1); -g*t; zeros(5, 1)];
    else
        x_set = [a*exp(-b*t); zeros(5, 1); -a*b*exp(-b*t); zeros(5, 1)];
    end
end