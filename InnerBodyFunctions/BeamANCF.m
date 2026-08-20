function [K_loc,Fe] = BeamANCF(Body,k,CalculateStiffness,FiniteDiference)
     
    xloc = Body.xloc;
    ElementDofs = Body.ElementDofs;
    
    K_loc = zeros(ElementDofs,ElementDofs);  % initialization of local stiffness matrix
    uk=Body.u(xloc(k,:));              
    qk0=Body.q0(xloc(k,:)); 

    Gint = Body.Gint;
    Nint = Body.Nint;    
    
    switch FiniteDiference

           case {"Matlab", "Matlab_automatic"}
               
                sqrtEps = sqrt(eps);
                
                qk=Body.q(xloc(k,:));   
                Phik=Body.Phim(k,:)';
                
                S = @(F_,xi,eta,zeta) Body.S(F_,qk0,Phik,xi,eta,zeta);  
                F = @(q,u,xi,eta,zeta) Body.F(q,qk0,u,xi,eta,zeta);
                
                dEdq =  Body.dEdq;
                dF_dq = @(xi,eta,zeta) Body.dF_dq(qk0,xi,eta,zeta);

                Fe=Fe_fun(qk,uk,F,dF_dq,dEdq,S,ElementDofs,Gint,Nint);

                if CalculateStiffness

                    if FiniteDiference == "Matlab" 
    
                       h = 2*sqrtEps;
                       Feh_all = zeros(ElementDofs, ElementDofs);
                       H = diag(h * ones(1, ElementDofs));
                       for jj = 1:ElementDofs
                            ukh = uk - H(:,jj);
                            qkh = qk - H(:,jj);   
                            Feh_all(:,jj) = Fe_fun(qkh,ukh,F,dF_dq,dEdq,S,ElementDofs,Gint,Nint);
                       end
        
                       K_loc = (Fe - Feh_all) ./ diag(H)';
                         
                    else % Body.FiniteDiference == "Matlab_automatic" 
    
                       fac = 1e-6;
                       if Body.SolutionBase == "Position"                    
                           G = @(t,y) Fe_fun(y,uk,F,dF_dq,dEdq,S,ElementDofs,Gint,Nint); 
                           K_loc = numjac(G, 0, qk, Fe, fac, []);
                       else
                           G = @(t,y) Fe_fun(qk,y,F,dF_dq,dEdq,S,ElementDofs,Gint,Nint);
                           K_loc = numjac(G, 0, uk, Fe, fac, []);
                       end 
                    end
                end
                                
           case "AceGen"
                Dvec = Body.Dvec;                
                Fe0 = K_loc(:,1);
                DIM = Body.DIM;
                DofsAtNode = Body.DofsAtNode;
                qk0f=Body.q0f(xloc(k,:));
                % Reshaping to adjust for AceGen
                qk0f_DIM = reshape(qk0f, [DIM, DofsAtNode])';
                qk0_DIM = reshape(qk0, [DIM, DofsAtNode])';
                uk_DIM = reshape(uk, [DIM, DofsAtNode])'; 
                [~,~,~,~,K_loc,Fe,~,~] = AceGenForce(qk0f_DIM,qk0_DIM,uk_DIM,Dvec,K_loc,Fe0,Gint',Nint);
                Fe = - Fe; % taking into account the difference between AceGen and Fe_fun                                  
    end           
