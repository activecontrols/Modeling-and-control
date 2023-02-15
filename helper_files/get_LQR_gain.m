function K = get_LQR_gain(A, B, Q, R)
    [K, ~, ~] = lqr(A, B, Q, R);
    Z = [A -(B/R*B'); -Q -A'];
    [U, S] = schur(Z);
    [U, S] = ordschur(U, S, 'lhp');
    [m,n] = size(U);
    U11 = U(1:(m/2), 1:(n/2));
    U21 = U((m/2+1):m, 1:(n/2));
    P = U21/U11;
    Kcheck = inv(R)*B'*P;
    disp(Kcheck - K);
end