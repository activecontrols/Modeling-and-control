close all
clear
clc

data = readmatrix("C:\Users\jamor\OneDrive\Desktop\imu.csv");

omega = data(:,6);  %rad/sec
time = (data(:,1) - data(1,1))/1000;  %sec

Fs = 1 / mean(diff(time));

save('LoadSingleAxisGyroscope', 'omega', 'Fs');