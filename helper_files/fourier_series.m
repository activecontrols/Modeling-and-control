function [F_t, L] = fourier_series(coarse_t, coarse_signal, num_coefficients, num_points)

time = linspace(coarse_t(1), coarse_t(end), num_points);
signal = interp1(coarse_t, coarse_signal, time);

L = time(end);
numCoeff = num_coefficients;

A0 = (1/2/L) * trapz(time, signal);

An = zeros(1, numCoeff);
Bn = zeros(1, numCoeff);

syms t F_t;
F_t(t) = A0;

for n = 1:numCoeff
    weightedCos = zeros(1, length(time));
    weightedSin = zeros(1, length(time));
    for m = 1:length(time)
        weightedCos(m) = signal(m) * cos(n * pi * time(m) / L);
        weightedSin(m) = signal(m) * sin(n * pi * time(m) / L);
    end
    An(n) = trapz(time, weightedCos) / L;
    Bn(n) = trapz(time, weightedSin) / L;
    F_t = F_t + An(n)*cos(n*pi*t/L) + Bn(n)*sin(n*pi*t/L);
end

end