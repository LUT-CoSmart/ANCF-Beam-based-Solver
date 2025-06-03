function Body = AddTensors(Body)
% Inner Force functions
addpath('InnerForceFunctions');

% Add shape function
addpath('TensorDerivations\' + Body.ElementType + '\Matlab\ShapeFunctions');
path = 'TensorDerivations\' + Body.ElementType + '\' + Body.FiniteDiference;
switch Body.FiniteDiference 
   case "AceGen"   
        disp("For chosen finite difference scheme, deformations are ony finite and displacement-based") 
        path = path + '\' + Body.ElementName;
        if ~isfolder(path)
           error('This element is not yet implemnetewd in AceGen');
        end
        Body.DeformationType = "Finite";
        Body.SolutionBase = "Displacement";

   case "Matlab"   
        addpath('MaterialsSecondKirhhoff'); % functions to calculate the 2nd Kirhhoff tensors in Matlab

        if Body.Fibers
           addpath(path + '\FiberVectors');
        end

        path = path + '\' + Body.SolutionBase;
        if Body.SolutionBase == "Position"
           disp("For chosen finite difference and solution-based scheme, deformations are only finite")
           Body.DeformationType = "Finite";
        else    
           path = path + '\' + Body.DeformationType; 
        end
        
   otherwise
        error('****** Choose correct Finite Diference scheme ******\n')
end        
 

addpath(path);