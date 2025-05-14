function xlocAll = xlocAllBeam(Body)

% function xlocAll makes full xloc for all elements 
% row - element, column - dof
xlocAll = zeros(Body.ElementNumber,Body.ElementDofs);

nloc = Body.nloc;
DofsAtNode = Body.DofsAtNode;
for k = 1:Body.ElementNumber
    loc = [];
    for j = 1: Body.ElementNodes        
        loc = [loc (nloc(k,j)-1)*DofsAtNode+1:nloc(k,j)*DofsAtNode];
    end    
     xlocAll(k,:) = loc;
end    