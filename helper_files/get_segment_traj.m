function [xsegment, usegment, tsegment, Ksegment, root] = get_segment_traj(segArray)
MOI = segArray{1};
constants = [segArray{3}; segArray{2}; segArray{4}];
xcrit1 = segArray{5}(:, 1);
xcrit2 = segArray{5}(:, 2);
ucrit2 = segArray{6};
limits = segArray{7};
throttleConsts = segArray{9};
numPoints = segArray{10};
ti = segArray{11};
tf = segArray{12};
Qbry = segArray{13};
Rbry = segArray{14};
genOn = segArray{15};

%Get symbolic EOMS and variables
[x, u, ~, Jx, Ju, consts, Jmat] = EOMS(throttleConsts);
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

%Determine lqr cost functions Q, and R
if genOn
    [Q, R, root] = genetic_algorithm(x, A, B, segArray);
else
    Q = diag(Qbry);
    R = diag(Rbry);
end

%Optimal control gain matrix K, solution S, and poles P
try
    [Ksegment, ~, ~] = lqr(A, B, Q, R);
catch e
    disp(e.message);
    error("LQR gain generation threw the error above!");
end
    
%Simulate using input data
%Utilizes dynamics' translational symmetry to approach critical points
[xsegment, usegment, tsegment] = simulate(ti, tf, numPoints, Ksegment, constants, MOI, xcrit1-xcrit2, limits);

% Create tracking gains for the simulation
C = [eye(3), zeros(3, length(A)-3)];
Atilde = [A, zeros(size(A, 1), 3); C, zeros(size(C, 1), 3)];
Btilde = [B; zeros(3, size(B, 2))];
Qtilde = diag([diag(Q); (.01*ones(3, 1)).^-2]);
Rtilde = R; 

[Ksegment, ~, ~] = lqr(Atilde, Btilde, Qtilde, Rtilde);

end