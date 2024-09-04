%clc, clear all, close all;
divv=pi/4;
divv2=pi/6;
Phi1=pi/3;
R1=1;
D1=2*R1;
%% 1
phi=pi:divv:3*pi;
c_x1=D1*cos(0*Phi1);c_y1=D1*sin(0*Phi1);
curve1a=[c_x1+R1*cos(phi)]';curve1b=[c_y1+R1*sin(phi)]';curve1=[curve1a curve1b];
%% 2
phi2=0:divv2:pi/3;
curve2a=[R1*cos(phi2)]';curve2b=[R1*sin(phi2)]';curve2=[curve2a curve2b];
%% 3
phi=4*pi/3:divv:2*pi+4*pi/3;
c_x3=D1*cos(1*Phi1);c_y3=D1*sin(1*Phi1);
curve3a=[c_x3+R1*cos(phi)]';curve3b=[c_y3+R1*sin(phi)]';curve3=[curve3a curve3b]; 
%% 4
phi2=pi/3:divv2:2*pi/3;
curve4a=[R1*cos(phi2)]';curve4b=[R1*sin(phi2)]';curve4=[curve4a curve4b];
%% 5
phi=5*pi/3:divv:2*pi+5*pi/3;
c_x5=D1*cos(2*Phi1);c_y5=D1*sin(2*Phi1);
curve5a=[c_x5+R1*cos(phi)]';curve5b=[c_y5+R1*sin(phi)]';curve5=[curve5a curve5b]; 
%% 6
phi2=2*pi/3:divv2:pi;
curve6a=[R1*cos(phi2)]';curve6b=[R1*sin(phi2)]';curve6=[curve6a curve6b];
%% 7
phi=0:divv:2*pi;
c_x7=D1*cos(3*Phi1);c_y7=D1*sin(3*Phi1);
curve7a=[c_x7+R1*cos(phi)]';curve7b=[c_y7+R1*sin(phi)]';curve7=[curve7a curve7b];
%% 8
phi2=pi:divv2:4*pi/3;
curve8a=[R1*cos(phi2)]';curve8b=[R1*sin(phi2)]';curve8=[curve8a curve8b];
%% 9
phi=pi/3:divv:2*pi+pi/3;
c_x9=D1*cos(4*Phi1);c_y9=D1*sin(4*Phi1);
curve9a=[c_x9+R1*cos(phi)]';curve9b=[c_y9+R1*sin(phi)]';curve9=[curve9a curve9b];
%% 10
phi2=4*pi/3:divv2:5*pi/3;
curve10a=[R1*cos(phi2)]';curve10b=[R1*sin(phi2)]';curve10=[curve10a curve10b];
%% 11
phi=2*pi/3:divv:2*pi+2*pi/3;
c_x11=D1*cos(5*Phi1);c_y11=D1*sin(5*Phi1);
curve11a=[c_x11+R1*cos(phi)]';curve11b=[c_y11+R1*sin(phi)]';curve11=[curve11a curve11b];
%% 12
phi2=5*pi/3:divv2:2*pi;
curve12a=[R1*cos(phi2)]';curve12b=[R1*sin(phi2)]';curve12=[curve12a curve12b];
% collect all original data
data_1=[curve1;curve2;curve3;curve4;curve5;curve6;curve7;curve8;curve9;curve10;curve11;curve12];
% plot(data_1(:,1),data_1(:,2),'.k');
% %% To be sure that all points in the same direction (counterclockwise)
% hold on
% for n = 1:length(data_1(:,1))
%     text(data_1(n,1),data_1(n,2),num2str(n))   
% end
%%
max_x=max(data_1(:,1));
min_x=min(data_1(:,1));
max_y=max(data_1(:,2));
min_y=min(data_1(:,2));

curve1a=[change(curve1(:,1),min_x,max_x) change(curve1(:,2),min_y,max_y)];
curve2a=[change(curve2(:,1),min_x,max_x) change(curve2(:,2),min_y,max_y)];
curve3a=[change(curve3(:,1),min_x,max_x) change(curve3(:,2),min_y,max_y)];
curve4a=[change(curve4(:,1),min_x,max_x) change(curve4(:,2),min_y,max_y)];
curve5a=[change(curve5(:,1),min_x,max_x) change(curve5(:,2),min_y,max_y)];
curve6a=[change(curve6(:,1),min_x,max_x) change(curve6(:,2),min_y,max_y)];
curve7a=[change(curve7(:,1),min_x,max_x) change(curve7(:,2),min_y,max_y)];
curve8a=[change(curve8(:,1),min_x,max_x) change(curve8(:,2),min_y,max_y)];
curve9a=[change(curve9(:,1),min_x,max_x) change(curve9(:,2),min_y,max_y)];
curve10a=[change(curve10(:,1),min_x,max_x) change(curve10(:,2),min_y,max_y)];
curve11a=[change(curve11(:,1),min_x,max_x) change(curve11(:,2),min_y,max_y)];
curve12a=[change(curve12(:,1),min_x,max_x) change(curve12(:,2),min_y,max_y)];
% collect changed data
data_2=[curve1a;curve2a;curve3a;curve4a;curve5a;curve6a;curve7a;curve8a;curve9a;curve10a;curve11a;curve12a];
data={curve1a,curve2a,curve3a,curve4a,curve5a,curve6a,curve7a,curve8a,curve9a,curve10a,curve11a,curve12a};
nu2=max(size(data));