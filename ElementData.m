% Define element properties
if Element == 3243
    Slope_x = true;  % Does the element include slope/derivative along x axis  
elseif Element == 3333 
    Slope_x = false;
elseif Element == 3343
    Slope_x = false;
elseif Element == 3353
    Slope_x = false;
elseif Element == 3363
    Slope_x = false;
elseif Element == 34103
    Slope_x = false;
else    
    fprint('Choose write element type!!!!')
end 
ElementName = num2str(Element);               % using 'abcd' classification, see in https://doi.org/10.1007/s11071-022-07518-z
ElemNodes = str2double(ElementName(2));       % Number of nodes            
DIM = str2double(ElementName(end));           % Problem dimensionality     
VecAtNode = str2double(ElementName(3:end-1)); % Vector functions per node
DofsAtNode = DIM * VecAtNode;                 % Number of Dofs in each element node
ElemDofs=ElemNodes*DofsAtNode;                % Total number of Dofs in the chosen element
