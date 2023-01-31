function traj = calc_trajectory(timeVector, xInit, linAcc, angAcc)
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[timeVector, x] = ode45(@(t,x) fun(t, x, linAcc, angAcc), [timeVector(1), timeVector(end)], xInit, options);
fplot(angAcc(2), [timeVector(1) timeVector(end)]);
figure;
plot(timeVector, x(:, 5));

traj = x;

end

function dxdt = fun(t, x, linAcc, angAcc)
    dxdt(1:6) = x(7:12)';
    linAccEval = subs(linAcc, t);
    angAccEval = subs(angAcc, t);
    dxdt(7:12) = [linAccEval; angAccEval];
    dxdt = dxdt';
end