%%  1 "edge" circular data 
%% S=0.075; 
R1=0.1;
angle=0:15:180;
curve1=R1*[ cosd(angle);sind(angle)]';
angle=180:15:360;
curve2=R1*[cosd(angle);sind(angle)]';    
% collect all original data
data_1=[curve1;curve2];
%%
max_x=max(data_1(:,1));
min_x=min(data_1(:,1));
max_y=max(data_1(:,2));
min_y=min(data_1(:,2));
% 
curve1a=[change(curve1(:,1),min_x,max_x) change(curve1(:,2),min_y,max_y)];
curve2a=[change(curve2(:,1),min_x,max_x) change(curve2(:,2),min_y,max_y)];
% collect changed data
data_2=[curve1a;curve2a];
data={curve1a,curve2a};
nu2=max(size(data));