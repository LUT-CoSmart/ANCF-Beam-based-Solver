% Loading conditions 
if Case == 1 % elongation load
   Fx_max = 9500;   
   Fx=Fx_max/steps:Fx_max/steps:Fx_max; 
   Fy=zeros(1,length(Fx));
   Fz=zeros(1,length(Fx));
   BC=0;       % 0 - reduced, 1 - full 
elseif Case == 2 % bending load
   Fy_max = 62.5*10^3; %5e5*0.5^8;   
   Fy=Fy_max/steps:Fy_max/steps:Fy_max; 
   Fx=zeros(1,length(Fy));
   Fz=zeros(1,length(Fy));
   BC=1;       % 0 - reduced, 1 - full 
else 
   disp('****** The loading type is not recognized ******');
   return;
end    
Fmesh=[Fx; Fy; Fz];