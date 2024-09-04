% Geometrical setup
addpath('CrossSections');
L=1; % Length
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

