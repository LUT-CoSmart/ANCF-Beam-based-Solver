%% Data 
R1=1;
R2=R1/2;

n = 10;
% Curve 1: R1, from -90° to 90° (bottom-left to top-right arc)
theta1 = linspace(-pi/2, pi/2, n);
curve1 = [R1*cos(theta1)', R1*sin(theta1)'];

% Curve 2: straight radial inward from R1 to R2 at 90°
r2 = linspace(R1, R2, n/2);
curve2 = [r2'.*cos(pi/2), r2'.*sin(pi/2)];

% Curve 3: R2, from 90° to -90° (top-right to bottom-right arc)
theta3 = linspace(pi/2, -pi/2, n);
curve3 = [R2*cos(theta3)', R2*sin(theta3)'];

% Curve 4: straight radial outward from R2 to R1 at -90°
r4 = linspace(R2, R1, n/2);
curve4 = [r4'.*cos(-pi/2), r4'.*sin(-pi/2)];

% collect all original data
data_1={curve1;curve2;curve3;curve4};
