function [x_set, t1] = get_trajectory(t)
    h0 = 60;
    hf = 1;
    g = 9.81;
    t1 = 1;
    b = (g*t1)/(h0 - hf - 0.5*g*t1^2);
    a = g*t1/(b*exp(-b*t1));
    if (t < t1)
        x_set = [h0 - (g/2)*t^2; zeros(5, 1); -g*t; zeros(5, 1)];
    else
        x_set = [a*exp(-b*t) + hf; zeros(5, 1); -a*b*exp(-b*t); zeros(5, 1)];
    end
end