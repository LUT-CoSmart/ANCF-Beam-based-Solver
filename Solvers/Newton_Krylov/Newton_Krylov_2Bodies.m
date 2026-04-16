function [u_bc,deltaf,Gap] = Newton_Krylov_2Bodies(ii, Body1,Body2,ContactType,ContactVariable,ContactFiniteDiference,Fext,Re,name)
        
         if nargin < 9
            name = "JF";
            warning('chosen JF algorithm');  
         end
         
         persistent v r s Fe Fc J;

         bc = [Body1.bc Body2.bc];

         q1_backup = Body1.q; % restore
         q2_backup = Body2.q; % restore
         
         u1_backup = Body1.u; % restore
         u2_backup = Body2.u; % restore
        
         q_backup = [q1_backup; q2_backup];

         Gap = 0;

         if ii == 1

            [Ke1,Fe1] = InnerForce(Body1); 
            [Ke2,Fe2] = InnerForce(Body2); % inner forces of the second body
            [Kc,Fc,Gap] = Contact(Body1,Body2,ContactType,ContactVariable,ContactFiniteDiference); % Contact forces
            
            Fe = [Fe1; Fe2];
            r =  -(Fe - Fext + Fc);
            s = zeros(size(r));
            v = r;
            u_bc =  zeros(size(r(bc)));

            % assembley
            Ke = [Ke1 zeros(Body1.TotalDofs,Body2.TotalDofs);
                  zeros(Body2.TotalDofs,Body1.TotalDofs) Ke2];            
            J = Kc + Ke;

         else
            u_bc =  zeros(size(r(bc)));
            nx = norm(q_backup(bc));        
            nv = norm(v(bc));
            h = 1e-10 * (1 + nx) / max(nv, 1e-20);

            r_bc = r(bc);
            v(~bc) = 0;
            v_bc = v(bc);

            switch name 
    
                case "JF"
                    g = J*v;
                    
                    Body1.SurfacePoints = feval(Body1.SurfacefunctionName, Body1, Body1.q);
                    Body2.SurfacePoints = feval(Body2.SurfacefunctionName, Body2, Body2.q);

                    Outcome = FindProjection(Body1.SurfacePoints, Body1.IsoData, Body2);

                    if ~isempty(Outcome)   
                        Gap = abs(sum(Outcome(:,5))); 
                    end

                    

                case "CG"
                    % Separation
                    Body1.u = Body1.u + h*v(1:Body1.TotalDofs);
                    Body1.q = Body1.q + h*v(1:Body1.TotalDofs);
            
                    Body2.u = Body2.u + h*v(Body1.TotalDofs + 1:end);        
                    Body2.q = Body2.q + h*v(Body1.TotalDofs + 1:end);

                    [~,Fe1] = InnerForce(Body1); 
                    [~,Fe2] = InnerForce(Body2); % inner forces of the second body
                    Fev = [Fe1; Fe2];

                    CalculateStiffness = false;
                    [~,Fcv,Gap] = Contact(Body1,Body2,ContactType,ContactVariable,ContactFiniteDiference,CalculateStiffness); % Contact forces
                    
                    
                    g = (Fev - Fe + Fcv - Fc) / h;

                    Body1.q = q1_backup; % restore
                    Body2.q = q2_backup; % restore
    
                    Body1.u = u1_backup; % restore
                    Body2.u = u2_backup; % restore

                otherwise
                    warning(' not correct Newton-Krylov algorithm, switched to JF');    
            end
    
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