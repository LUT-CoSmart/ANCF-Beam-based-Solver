% Define element properties
addpath('ElementTensors\Matlab','ElementTensors\AceGen');
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

    