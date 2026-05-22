function [u_bc,deltaf,Gap] = Newton_Krylov_2Bodies(ii, Body1,Body2,ContactType,ContactVariable,ContactFiniteDiference,Fext,Re,name)
        

         % In many cases bodies aren't interesced, so stiffness matrix of the contact part is zero,
         % in this regard JF modification for Newton-Krylov is not working properly, and therefore, removed.

         persistent v r s Fe Fc;

         bc = [Body1.bc Body2.bc];

         q1_backup = Body1.q; % save
         q2_backup = Body2.q; % save
         
         u1_backup = Body1.u; % save
         u2_backup = Body2.u; % save
        
         q_backup = [q1_backup; q2_backup];

         if ii == 1

            [~,Fe1] = InnerForce(Body1); 
            [~,Fe2] = InnerForce(Body2); % inner forces of the second body
            [~,Fc,Gap] = Contact(Body1,Body2,ContactType,ContactVariable,ContactFiniteDiference,false); % Contact forces
    
            Fe = [Fe1; Fe2];
            r =  -(Fe - Fext + Fc);
            s = zeros(size(r));
            v = r;
            u_bc =  zeros(size(r(bc)));

         else
            u_bc =  zeros(size(r(bc)));
            nx = norm(q_backup(bc));        
            nv = norm(v(bc));
            h = 1e-7 * (1 + nx) / max(nv, 1e-20);

            r_bc = r(bc);
            v(~bc) = 0;
            v_bc = v(bc);

            % Separation
            Body1.u = Body1.u + h*v(1:Body1.TotalDofs);
            Body1.q = Body1.q + h*v(1:Body1.TotalDofs);
            
            Body2.u = Body2.u + h*v(Body1.TotalDofs + 1:end);        
            Body2.q = Body2.q + h*v(Body1.TotalDofs + 1:end);

            [~,Fe1] = InnerForce(Body1); 
            [~,Fe2] = InnerForce(Body2); % inner forces of the second body
            Fev = [Fe1; Fe2];

            [~,Fcv,Gap] = Contact(Body1,Body2,ContactType,ContactVariable,ContactFiniteDiference,false); % Contact forces
                                        
            g = (Fev - Fe + Fcv - Fc) / h;

            Body1.q = q1_backup; % restore
            Body2.q = q2_backup; % restore
    
            Body1.u = u1_backup; % restore
            Body2.u = u2_backup; % restore
    
            g_bc = g(bc);
            
            alpha = r_bc' * r_bc/(v_bc'*g_bc);
            s = s + alpha * v;
            r = r - alpha * g; % update r
    
            beta = r(bc)' * r(bc) / (r_bc' *r_bc);
            v = r + beta * v;

        end

        deltaf=r(bc)/norm(Fext(bc));

        if all(abs(deltaf) < Re) % to make the exit from the algorith similar to others
           u_bc = s(bc);
        end   