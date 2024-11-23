% Material parameters
if (Material==1) || (Material==2) || (Material==5) || (Material==0)    
   if Material==1 % Neo-Hookean
      Fibers = false;  % Material is isotropic
      mu=9e5;
      const=mu;      
      MaterialName = 'Neo';      
   elseif Material==2 % 2 contant Mooney-Rivlin
      Fibers = false;
      c10=33.4e4;
      c01 = -337;
      const=[c10, c01];  
      MaterialName = 'Mooney2';
   elseif Material==5 % 5-contant Mooney-Rivlin 
      Fibers = false;
      c10 = -7.7e5;
      c01 = 9.1e5;
      c11 = 1.03e6;
      c20 = -2.7e5;
      c02 = -5.9e5;
      const=[c10, c01, c11, c20, c02];
      MaterialName = 'Mooney5';
   elseif Material==0 % GOH   
      Fibers = true;
      c10 = 7.64e3;
      k1 = 996.6e3;
      k2 = 524.6;
      a0 = [1 0 0];                % fiber direction
      kappa = 0;                   % fiber dipersion
      const=[c10, k1, k2, kappa, a0]; 
      MaterialName = 'GOH';      
   end
   d = 10^(-10); % bulk module
else
   if Material==3 % K.-S.
      Fibers = false;
      E=2.07e11;
      nu=0.3;
      const = [E,nu];
      MaterialName = 'KS';
      d =[];     
   else   
     disp('****** The material type is not recognized ******');
     return;
   end  
end
if Fibers
   fiber_twist = 0;  % inner (fiber) pre-twist 
else
   fiber_twist = 0;
end    
const = [const, d];