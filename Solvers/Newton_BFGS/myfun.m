function [ff_bc,K_bc] = myfun(Body,u_bc)

        Fext = Body.Fext;
        

        Body.u(Body.bc) = Body.u(Body.bc)+u_bc; % Add displacement to previous one
        Body.q(Body.bc) = Body.q(Body.bc)+u_bc; % change the global positions

        [K,Fe] = InnerForce(Body);
        ff =  Fe - Fext;                 % assembley
        
        K_bc = K(Body.bc,Body.bc);
        ff_bc= ff(Body.bc);              % Eliminate linear constraints from force vector        
        