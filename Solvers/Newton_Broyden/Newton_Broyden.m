function [u_bc,deltaf] = Newton_Broyden(step, Body, Fext)
       
       persistent Bn_m1 ff_bc_old u_bc_old
        
       [~,Fe] = InnerForce(Body);
       ff =  Fe - Fext;                 % assembley
       ff_bc=ff(Body.bc);               % Eliminate linear constraints from force vector 
        
       if step == 1  
            [K,~] = InnerForce(Body);        
            Bn_m1 = inv( K(Body.bc,Body.bc) );
       else   
            zn = ff_bc - ff_bc_old;
            s = u_bc_old;

            denom = s' * Bn_m1 * zn;
            Bn_m1 = Bn_m1 + (s - Bn_m1 * zn ) * s' * Bn_m1 ./ denom ;
       end       

       u_bc = - Bn_m1*ff_bc;
        
       deltaf=ff_bc/norm(Fext(Body.bc));

       % Saving for the next iteration
       ff_bc_old = ff_bc;
       u_bc_old = u_bc;
