% Material parameters
addpath('MaterialsSecondKirhhoff');
if (Material==1) || (Material==2) || (Material==5) || (Material==0)    
   if Material==1 % Neo-Hookean
      mu=9e5;
      const=mu;
      MaterialName = 'Neo';
   elseif Material==2 % 2 contant Mooney-Rivlin
      c10=33.4e4;
      c01 = -337;
      const=[c10, c01];  
      MaterialName = 'Mooney2';
   elseif Material==5 % 5-contant Mooney-Rivlin   
      c10 = -7.7e5;
      c01 = 9.1e5;
      c11 = 1.03e6;
      c20 = -2.7e5;
      c02 = -5.9e5;
      const=[c10, c01, c11, c20, c02];
      MaterialName = 'Mooney5';
   elseif Material==0 % GOH   
      c10 = 7.64e3;
      k1 = 996.6e3;
      k2 = 524.6;
      kappa = 0;
      a = [sqrt(1)/2 0 sqrt(3)/2];
      const=[c10, k1, k2, kappa, a]; 
      MaterialName = 'GOH';
   end
   d = 10^(-10); % bulk module for incompresible materials
   const = [const, d];
else
   if Material==3 % K.-S.
      E=2.07e11;
      nu=0.3;
      const = [E,nu];
      MaterialName = 'KS';
   else   
     disp('****** The material type is not recognized ******');
     return;
   end  
end
    