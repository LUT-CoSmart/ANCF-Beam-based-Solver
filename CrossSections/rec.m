%% Rectangular
hh=0.05;
x=-0.05:hh:0.05;
y=-0.05*ones(1,length(x));
curve1=[x' y'];
y=-0.05:hh:0.05;
x=0.05*ones(1,length(y));
curve2=[x' y'];
x=0.05:-hh:-0.05;
y=0.05*ones(1,length(x));
curve3=[x' y'];
y=0.05:-hh:-0.05;
x=-0.05*ones(1,length(y));
curve4=[x' y'];

% collect all original data
data_1=[curve1;curve2;curve3;curve4];

max_x=max(data_1(:,1));
min_x=min(data_1(:,1));
max_y=max(data_1(:,2));
min_y=min(data_1(:,2));

curve1a=[change(curve1(:,1),min_x,max_x) change(curve1(:,2),min_y,max_y)];
curve2a=[change(curve2(:,1),min_x,max_x) change(curve2(:,2),min_y,max_y)];
curve3a=[change(curve3(:,1),min_x,max_x) change(curve3(:,2),min_y,max_y)];
curve4a=[change(curve4(:,1),min_x,max_x) change(curve4(:,2),min_y,max_y)];
% collect changed data
data_2=[curve1a;curve2a;curve3a;curve4a];
data={curve1a,curve2a,curve3a,curve4a};
nu2=max(size(data));