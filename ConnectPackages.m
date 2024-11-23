%% Connect necessary packages
% Support functions
addpath('SupportFunctions');
% For Geometry approximation
addpath('MeshFunctions')
addpath('GaussPoints');
%% Problem data
addpath('ProblemData'); 
% Files include data for all problems (required to be activated)
Materials;
Loads;
Geometry;
Approximation;
%% Connect to inner energy functions
% this document allows, based on the chosen options, choose the right folder 
% for computing inner energy.
ElementName = string(Element);
% Inner Force calculation
addpath('InnerForceFunctions');
if sol_acegen
   path = 'TensorDerivations/AceGen/' + ElementName; 
   addpath(path);
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
   if Fibers
      addpath('TensorDerivations/Matlab/FiberVectors');
   end    
end  
