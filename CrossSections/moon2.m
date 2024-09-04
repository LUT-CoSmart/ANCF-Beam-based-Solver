%%  1 "edge" circular data 
%% S=0.075; 
R1=0.233596680327607;
R2=R1/2;
curve1=[    R1*cos(pi/2)      -R1*sin(pi/2);
            R1*cos(pi/4)      -R1*sin(pi/4);
               R1*cos(0)          R1*sin(0);
            R1*cos(pi/4)       R1*sin(pi/4);
            R1*cos(pi/2)       R1*sin(pi/2);];
        
curve2=[    R1*cos(pi/2)      R1*sin(pi/2);
       (R1+R2)/2*cos(pi/2) (R1+R2)/2*sin(pi/2);
            R2*cos(pi/2)      R2*sin(pi/2);];
        
        
curve3=[    R2*cos(pi/2)      R2*sin(pi/2);
            R2*cos(pi/4)      R2*sin(pi/4);
               R2*cos(0)         R2*sin(0);
            R2*cos(pi/4)     -R2*sin(pi/4);
            R2*cos(pi/2)     -R2*sin(pi/2);];
        
        
curve4=[     R2*cos(pi/2)      -R2*sin(pi/2);
       -(R1+R2)/2*cos(pi/2) -(R1+R2)/2*sin(pi/2);   
            -R1*cos(pi/2)      -R1*sin(pi/2);];
% collect all original data
data_1=[curve1;curve2;curve3;curve4];
%%
max_x=max(data_1(:,1));
min_x=min(data_1(:,1));
max_y=max(data_1(:,2));
min_y=min(data_1(:,2));
% 
S=(max_y/2-min_y/2)*(max_x/2-min_x/2);
W=max_x-min_x;
H=max_y-min_y;

curve1a=[change(curve1(:,1),min_x,max_x) change(curve1(:,2),min_y,max_y)];
curve2a=[change(curve2(:,1),min_x,max_x) change(curve2(:,2),min_y,max_y)];
curve3a=[change(curve3(:,1),min_x,max_x) change(curve3(:,2),min_y,max_y)];
curve4a=[change(curve4(:,1),min_x,max_x) change(curve4(:,2),min_y,max_y)];
% collect changed data
data_2=[curve1a;curve2a;curve3a;curve4a];
%plot(data_2(:,1),data_2(:,2),'ok')
data={curve1a,curve2a,curve3a,curve4a};
nu2=max(size(data));