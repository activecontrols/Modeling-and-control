m = 4; %mass in Kg
R = 95; %Radius in mm
h = 1; %height in m

R = R * 0.001;
moi_tensor = [ [m/4*(1/3*h^2+R^2) 0 0]; [0 m/4*(1/3*h^2+R^2) 0]; [0 0 (m/2*R^2)]];