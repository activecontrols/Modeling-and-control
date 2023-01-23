numGraphs = 6;
fileName = 'Loki_L840CT.csv';
raw = readmatrix(fileName);
time = raw(:,1);
thrust = raw(:, 2);

for n=1:numGraphs
    subplot(2, 3, n);
    fplot(thrust_series(time, thrust, n^2), [time(1), time(end)]);
    plotTitle = sprintf("Series with %d coefficients", n^2);
    hold on;
    plot(time, thrust);
    title(plotTitle);
    xlabel("Time (s)");
    ylabel("Force (N)");
    legend("Fourier Series", "Real data");
end