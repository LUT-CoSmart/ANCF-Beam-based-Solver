function xlocAll= xlocAllANCF(ElemNodes,DofsAtNode,nloc)
% function xlocAllANCF makes full xloc for all elements 
% row - element, column - dof
[nl,~] = size(nloc);
xlocAll = zeros(nl,DofsAtNode*ElemNodes);
for k = 1:nl
    loc = [];
    for j = 1: ElemNodes        
        loc = [loc (nloc(k,j)-1)*DofsAtNode+1:nloc(k,j)*DofsAtNode];
    end    
     xlocAll(k,:) = loc;
end    