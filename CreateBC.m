% Creates boundary conditions
% create global vector of nodal coordinates
for jj=1:nn
   q((jj-1)*DofsAtNode+1:jj*DofsAtNode)=P0(jj,:); 
end
% Define vector of linear constraints
bc = logical(ones(1,nx));
if BC == 0    
    if Element==3243
        Dofs = [1:3,5:7,9:11];
    elseif Element==3333
        Dofs = [1:3,4,6,7,8];
    elseif Element==3343
        Dofs = [1:3,4,6,7,8,10:12];    
    elseif Element==3353
        Dofs = [1:3,4,6,7,8,10:15];
    elseif Element==3363
        Dofs = [1:3,4,6,7,8,10:18];     
    elseif Element==34103
        Dofs = [1:3,4,6,7,8,10:30]; 
    end
    %% TODO: Generalize the implementaion of BC == 0 w/o dependency on element name
    % Dofs = [1:DIM,2+DIM+length_rx+1,2*DIM+length_rx-1,]    
    % Dofs  =  1:DofsAtNode;
    % Dofs = [1:DIM, length_rx/DIM * (DIM+ (DIM-1:DIM)), DIM + length_rx/DIM * DIM + [1,3:DIM], DIM + length_rx/DIM * DIM + DIM + [1:DIM-1]];
    %     Dofs = Dofs(Dofs > 0)
    % 
    % % Dofs = 1:DIM;
    % % for i = 1:DIM % this one will be important for unification with Advanced media descriptions
    % %     Dofs = [Dofs,  ];
    % % end    
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

     