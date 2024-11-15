% Geometrical setup
addpath('CrossSections');
L=1; % Length
if App ~= 0 % If we use non-standard integratrion and need points 
    if Area==1 % if Approximation~=0 it will be substituted
        rec;
    elseif Area==2
        moon;
    elseif Area==0
        ten;
    elseif Area==3
        moon2;
    elseif Area==4
        moon3;     
    else
        disp('****** The area type is not recognized ******');
        return;
    end
    W=max_x-min_x;
    H=max_y-min_y;
else % if we use standard inegration, then the elemnet dimensions can be set
    W = 0.1;
    H = 0.1;
end    
% Configuration pre-twist
twist_angle = 0; % Twist around X axis in degrees
ro = 0;          % Distance till the rotatoin center