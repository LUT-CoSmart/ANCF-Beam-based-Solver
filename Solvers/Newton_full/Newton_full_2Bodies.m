function [u_bc,deltaf,Gap] = Newton_cont2Body_full(Body1,Body2,ContactType,ContactVariable,ContactFiniteDiference,Fext)
         
         bc = [Body1.bc Body2.bc];

         % Elastic forces
         [Ke1,Fe1] = InnerForce(Body1);    
         [Ke2,Fe2] = InnerForce(Body2);


         % Contact forces
         [Kc,Fc,Gap] = Contact(Body1,Body2,ContactType,ContactVariable,ContactFiniteDiference);
         
         % Assemblance
         Fe = [Fe1; Fe2];
         Ke = [Ke1 zeros(Body1.TotalDofs,Body2.TotalDofs);
               zeros(Body2.TotalDofs,Body1.TotalDofs) Ke2];   
         K = Kc + Ke;
         ff =  Fe - Fext + Fc; 

         % Calculations
         K_bc = K(bc,bc); 
         ff_bc = ff(bc);
         deltaf=ff_bc/norm(Fext(bc));         
         u_bc = Solving(K_bc,ff_bc); 

