function [u_bc,deltaf] = Newton_full(Body,Fext)
        
        [K,Fe] = InnerForce(Body);
  
        K_bc = K(Body.bc,Body.bc);            % Eliminate linear constraints from stiffness matrix
        ff =  Fe - Fext;

        ff_bc=ff(Body.bc);               % Eliminate linear constraints from force vector
        deltaf=ff_bc/norm(Fext(Body.bc));% Compute residua

        u_bc = Solving(K_bc,ff_bc);  
