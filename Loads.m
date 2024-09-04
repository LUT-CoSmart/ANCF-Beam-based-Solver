% Loading conditions 
if Case == 1 % elongation load
   Fx_max = 9500;
   
   Fx=Fx_max/steps:Fx_max/steps:Fx_max; 
   Fy=zeros(1,length(Fx));
   Fz=zeros(1,length(Fx));
elseif Case == 2 % bending load
   Fy_max = 5e5*0.5^3;
   
   Fy=Fy_max/steps:Fy_max/steps:Fy_max; 
   Fx=zeros(1,length(Fy));
   Fz=zeros(1,length(Fy));
else 
   disp('****** The loading type is not recognized ******');
   return;
end    
Fmesh=[Fx; Fy; Fz];