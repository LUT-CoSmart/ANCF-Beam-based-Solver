% Define element properties
if Element == 3243
    ElemType='3D two noded 24 dofs ANCF beam';  
    ElemNodes=2;
    DofsAtNode=12;
    ElemDofs=ElemNodes*DofsAtNode;
    DIM=3;
elseif Element == 3333
    ElemType='3D three noded 27 dofs ANCF beam';  
    ElemNodes=3;
    DofsAtNode=9;
    ElemDofs=ElemNodes*DofsAtNode;
    DIM=3;
elseif Element == 3353
    ElemType='3D three noded 45 dofs ANCF beam';  
    ElemNodes=3;
    DofsAtNode=15;
    ElemDofs=ElemNodes*DofsAtNode;
    DIM=3;     
elseif Element == 3363
    ElemType='3D three noded 54 dofs ANCF beam';  
    ElemNodes=3;
    DofsAtNode=18;
    ElemDofs=ElemNodes*DofsAtNode;
    DIM=3;  
elseif Element == 34103
    ElemType='3D four noded 120 dofs ANCF beam';  
    ElemNodes=4;
    DofsAtNode=30;
    ElemDofs=ElemNodes*DofsAtNode;
    DIM=3; 
else    
    fprint('Choose write element type!!!!')
end 
% Access to element tensors and force functions
if sol_acegen == true
   % addpath('ElementTensors\Matlab','ElementTensors\AceGen');
else
   if disp_based == false
      addpath('TensorDerivations\Matlab\ElementTensors\Position');  
   else
       if small_deformation == false
           addpath('TensorDerivations\Matlab\ElementTensors\Displacement\Finite');
       else
           addpath('TensorDerivations\Matlab\ElementTensors\Displacement\Small');
       end
   end
end   
PosDofs = []; % identefication of positional DoFs within the chosen element 
for i = 1:ElemNodes
    PosDofs = [PosDofs (i-1)*DofsAtNode+1:(i-1)*DofsAtNode+3];
end    
