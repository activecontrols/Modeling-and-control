function plotTrajectory(xPlot, uPlot, tScale, scaleFactor, numSkip, constants)
thrustArm = constants(2);
figure;
rTraj = xPlot(1:3, :);
rocketBody = zeros(6, length(tScale));
TVec = zeros(3, length(tScale));
maxThrust = max(uPlot(3, :));

for i = 1:length(tScale)
    rCurrent = rTraj(:, i);
    qshort = xPlot(7:9, i);
    q0C = sqrt(1 - qshort(1)^2 - qshort(2)^2 - qshort(3)^2);
    qCurrent = [q0C; qshort];
    CBICurrent = quatRot(qCurrent)';
    betaCurrent = uPlot(1, i);
    gammaCurrent = uPlot(2, i);
    F_TCurrent = .5*scaleFactor*uPlot(3, i)/maxThrust;
    TBCurrent = F_TCurrent * [cos(betaCurrent)*cos(gammaCurrent); cos(betaCurrent)*sin(gammaCurrent); -sin(betaCurrent)];
    rocketBody(1:3, i) = rCurrent + CBICurrent*[thrustArm*scaleFactor; 0; 0];
    rocketBody(4:6, i) = rCurrent + CBICurrent*[-thrustArm*scaleFactor; 0; 0];
    TVec(:, i) = rCurrent + CBICurrent*([-thrustArm*scaleFactor; 0; 0] - TBCurrent);
end
plot3(rTraj(2, :), rTraj(3, :), rTraj(1, :));
hold on;

for i = 1:numSkip:length(tScale)
    plot3([rocketBody(2, i); rocketBody(5, i)], [rocketBody(3, i); rocketBody(6, i)], [rocketBody(1, i); rocketBody(4, i)], 'b-');
    plot3([rocketBody(5, i); TVec(2, i)], [rocketBody(6, i); TVec(3, i)], [rocketBody(4, i); TVec(1, i)], 'r-');
end

title("Position Trajectory");
xlabel("R2 (m)");
ylabel("R3 (m)");
zlabel("R1 (m)");
axis square
grid on;

xMin = min([rocketBody(2, :); rocketBody(5, :); TVec(2, :)], [], 'all');
xMax = max([rocketBody(2, :); rocketBody(5, :); TVec(2, :)], [], 'all');
yMin = min([rocketBody(3, :); rocketBody(6, :); TVec(3, :)], [], 'all');
yMax = max([rocketBody(3, :); rocketBody(6, :); TVec(3, :)], [], 'all');
zMin = min([rocketBody(1, :); rocketBody(4, :); TVec(1, :)], [], 'all');
zMax = max([rocketBody(1, :); rocketBody(4, :); TVec(1, :)], [], 'all');

xRange = xMax - xMin;
yRange = yMax - yMin;
zRange = zMax - zMin;
xMid = (xMax + xMin) / 2;
yMid = (yMax + yMin) / 2;
zMid = (zMax + zMin) / 2;
plotRange = max([xRange, yRange, zRange]);

axis([xMid - 5*plotRange/8, xMid + 5*plotRange/8, yMid - 5*plotRange/8, yMid + 5*plotRange/8, zMid - 5*plotRange/8, zMid + 5*plotRange/8]);
end
function CIB = quatRot(q)
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    CIB = [1-2*(q2^2 + q3^2), 2*(q1*q2 + q0*q3), 2*(q1*q3 - q0*q2);
        2*(q1*q2 - q0*q3), 1-2*(q1^2 + q3^2), 2*(q2*q3 + q0*q1);
        2*(q1*q3 + q0*q2), 2*(q2*q3 - q0*q1), 1-2*(q1^2 + q2^2)];
end