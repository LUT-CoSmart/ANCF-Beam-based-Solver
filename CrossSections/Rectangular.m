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
data_1={curve1;curve2;curve3;curve4};