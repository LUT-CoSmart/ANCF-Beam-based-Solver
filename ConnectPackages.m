%% Connect packages
% Support functions
addpath('SupportFunctions');
% For Geometry approximation
addpath('MeshFunctions')
addpath('CrossSections');
addpath('GaussPoints');
% Inner Force calculation
addpath('InnerForceFunctions');
if sol_acegen
   addpath('TensorDerivations/AceGen');
else
   addpath('MaterialsSecondKirhhoff');  % functions to calculate the 2nd Kirhhoff tensot in Matlab
   if disp_based
      if small_deformation
         addpath('TensorDerivations\Matlab\ElementTensors\Displacement\Small');           
      else
         addpath('TensorDerivations\Matlab\ElementTensors\Displacement\Finite');
      end
   else
      addpath('TensorDerivations\Matlab\ElementTensors\Position');  
   end      
end  


