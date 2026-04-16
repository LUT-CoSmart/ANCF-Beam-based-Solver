function [u_bc,deltaf,Gap] = Newton_Broyden_2Bodies(ii, Body1,Body2,ContactType,ContactVariable,ContactFiniteDiference,Fext)
        
         bc = [Body1.bc Body2.bc];
         CalculateStiffness = false; 

         persistent Bn_m1 ff_bc_old u_bc_old
         
         % on each iteration we recalculate matrix
         [~,Fe1] = InnerForce(Body1); 
         [~,Fe2] = InnerForce(Body2); % inner forces of the second body
         [~,Fc,Gap] = Contact(Body1,Body2,ContactType,ContactVariable,ContactFiniteDiference,CalculateStiffness); % Contact forces
         
         Fe = [Fe1; Fe2];
         ff =  Fe - Fext + Fc;
         ff_bc = ff(bc);

         if ii == 1

            [Ke1,~] = InnerForce(Body1); 
            [Ke2,~] = InnerForce(Body2); % inner forces of the second body
            [Kc,~,~] = Contact(Body1,Body2,ContactType,ContactVariable,ContactFiniteDiference); % Contact forces
            
            % assembley
            Ke = [Ke1 zeros(Body1.TotalDofs,Body2.TotalDofs);
                  zeros(Body2.TotalDofs,Body1.TotalDofs) Ke2];            
            K = Kc + Ke;

            K_bc = K(bc,bc);  % applying dofs            
            Bn_m1 = inv(K_bc);

         else
            zn = ff_bc - ff_bc_old;
            s = u_bc_old;
            denom = s' * Bn_m1 * zn;
            Bn_m1 = Bn_m1 + (s - Bn_m1 * zn ) * s' * Bn_m1 ./ denom ;

         end

        u_bc = - Bn_m1*ff_bc;

        deltaf=ff_bc/norm(Fext(bc)); 

        % Saving for the next iteration
        ff_bc_old = ff_bc;
        u_bc_old = u_bc;
