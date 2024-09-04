function Fext=Fext_ANCF(DofsAtNode,nn,nx,F)
% Define vector of external forces
Fext = zeros(nx,1);
Fext(xlocANCF(DofsAtNode,nn,1)) = F(1);
Fext(xlocANCF(DofsAtNode,nn,2)) = F(2);
Fext(xlocANCF(DofsAtNode,nn,3)) = F(3);