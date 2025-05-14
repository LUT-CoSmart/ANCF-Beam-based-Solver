function Fext=Fext_ANCF(Body,F)
% Define vector of external forces
Fext = zeros(Body.nx,1);
Fext(xlocANCF(Body.DofsAtNode,Body.nn,1)) = F(1);
Fext(xlocANCF(Body.DofsAtNode,Body.nn,2)) = F(2);
Fext(xlocANCF(Body.DofsAtNode,Body.nn,3)) = F(3);