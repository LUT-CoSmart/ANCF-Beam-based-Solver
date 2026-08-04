function [K_loc,Fe] = BeamANCF(Body,k, CalculateStiffness)
     
    xloc = Body.xloc;
    ElementDofs = Body.ElementDofs;
    
    K_loc = zeros(ElementDofs,ElementDofs);  % initialization of local stiffness matrix
    uk=Body.u(xloc(k,:));              
    qk0=Body.q0(xloc(k,:)); 

    Gint = Body.Gint;
    Nint = Body.Nint;    
    switch Body.FiniteDiference

           case {"Matlab", "Matlab_automatic"}
               
                sqrtEps = sqrt(eps);
                DT = Body.DeformationType;
                SB = Body.SolutionBase;

                F0 = @(xi,eta,zeta) Body.F0(qk0,xi,eta,zeta);
                Nm_xi = Body.Nm_xi_;  
                                
                Fibers = Body.Fibers;
                ElemDim = 0.5 * Body.Length.Ln * Body.detF0;

                qk=Body.q(xloc(k,:));   
                Phik=Body.Phim(k,:)';
                
                if Fibers
                    a0_axis =  @(xi,eta,zeta) Body.a0_fib(qk0,Phik,xi,eta,zeta);
                    S = @(F_,xi,eta,zeta) Body.S(F_, a0_axis(xi,eta,zeta));  
                else
                    S = @(F_,xi,eta,zeta) Body.S(F_);  
                end

                F = @(q,u,xi,eta,zeta) Body.F(q,qk0,u,xi,eta,zeta);

                Fe=Fe_fun(Nm_xi,qk,uk,ElemDim,F,F0,S,ElementDofs,Gint,Nint,DT,SB);
                if CalculateStiffness

                    if Body.FiniteDiference == "Matlab" 
    
                       h = 2*sqrtEps;
                       Feh_all = zeros(ElementDofs, ElementDofs);
                       H = diag(h * ones(1, ElementDofs));
                       for jj = 1:ElementDofs
                            ukh = uk - H(:,jj);
                            qkh = qk - H(:,jj);   
                            Feh_all(:,jj) = Fe_fun(Nm_xi,qkh,ukh,ElemDim,F,F0,S,ElementDofs,Gint,Nint,DT,SB);
                       end
        
                       K_loc = (Fe - Feh_all) ./ diag(H)';
                         
                    else % Body.FiniteDiference == "Matlab_automatic" 
    
                       fac = 5e-7;
                       if SB == "Position"                    
                           G = @(t,y) Fe_fun(Nm_xi,y,uk,ElemDim,F,F0,S,ElementDofs,Gint,Nint,DT,SB); 
                           K_loc = numjac(G, 0, qk, Fe, fac, []);
                       else
                           G = @(t,y) Fe_fun(Nm_xi,qk,y,ElemDim,F,F0,S,ElementDofs,Gint,Nint,DT,SB);
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
