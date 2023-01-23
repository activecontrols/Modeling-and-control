function [F_t, L] = thrust_series(coarse_t, coarse_thrust, num_coefficients)

time = linspace(coarse_t(1), coarse_t(end), 1000);
thrust = interp1(coarse_t, coarse_thrust, time);

L = time(end);
numCoeff = num_coefficients;

A0 = (1/2/L) * trapz(time, thrust);

An = zeros(1, numCoeff);
Bn = zeros(1, numCoeff);

syms t F_t;
F_t(t) = A0;

for n = 1:numCoeff
    weightedCos = zeros(1, length(time));
    weightedSin = zeros(1, length(time));
    for m = 1:length(time)
        weightedCos(m) = thrust(m) * cos(n * pi * time(m) / L);
        weightedSin(m) = thrust(m) * sin(n * pi * time(m) / L);
    end
    An(n) = trapz(time, weightedCos) / L;
    Bn(n) = trapz(time, weightedSin) / L;
    F_t = F_t + An(n)*cos(n*pi*t/L) + Bn(n)*sin(n*pi*t/L);
end

end