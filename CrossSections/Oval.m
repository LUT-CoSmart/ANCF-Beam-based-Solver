
angle=0:15:360;
curve1=[ cosd(angle);sind(angle)]';
angle=180:15:360;
curve2=[cosd(angle);sind(angle)]';    
% collect all original data
data_1={curve1;curve2};