close all
clear
clc

data = readmatrix("C:\Users\jamor\OneDrive\Desktop\imu.csv");

time = (data(:,1) - data(1,1))/1000;  %sec
omega_v = data(:,3:5);  %rad/sec

Fs = 1 / mean(diff(time));

save('LoadSingleAxisGyroscope', 'omega_v', 'Fs');