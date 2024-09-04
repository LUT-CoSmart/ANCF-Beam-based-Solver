% Creates boundary conditions
% create global vector of nodal coordinates
for jj=1:nn
   ee((jj-1)*DofsAtNode+1:jj*DofsAtNode)=P0(jj,:); 
end
% Define vector of linear constraints
bc = logical(ones(1,nx));
if BC == 0    
    if Element==3243
        Dofs = [1:3,5:7,9:11];
    elseif Element==3333
        Dofs = [1:3,4,6,7,8];
    elseif Element==3353
        Dofs = [1:3,4,6,7,8,10:15];
    elseif Element==3363
        Dofs = [1:3,4,6,7,8,10:18];     
    elseif Element==34103
        Dofs = [1:3,4,6,7,8,10:30]; 
    end
elseif BC == 1
    Dofs = 1:DofsAtNode; 
else
    disp('****** The boundary conditions are not recognized ******');
    return;
end    
bcInd=xlocANCF(DofsAtNode,1,Dofs);
%bcInd
if bcInd~=0
    bc(bcInd)=0;% number of degrees of freedom of system after linear constraints  
end

     