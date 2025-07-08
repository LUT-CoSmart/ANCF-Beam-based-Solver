function [K_loc,Fe] = BeamANCF(Body,k)
     
    sqrtEps = sqrt(eps);
    xloc = Body.xloc;
    ElementDofs = Body.ElementDofs;
    K_loc = zeros(ElementDofs,ElementDofs);  % initialization of local stiffness matrix
    Fe0 = K_loc(:,1);
    Fe = Fe0;

    uk=Body.u(xloc(k,:));              
    qk0=Body.q0(xloc(k,:)); 

    Dvec = Body.Dvec;
    Gint = Body.Gint;
    Nint = Body.Nint;    

    switch Body.FiniteDiference
            
           case "Matlab"
                 h = 2*sqrtEps; 
                 
                 Fibers = Body.Fibers;
                 detF0 = Body.detF0;
                 PosDofs = Body.PosDofs;
                 qk=Body.q(xloc(k,:));
                 phik=Body.phim(k,:)';    
                 Phik=Body.Phim(k,:)';                 
                 Fe=Fe_fun(Fe0,uk,qk,qk0,phik,Phik,Fibers,Dvec,ElementDofs,PosDofs,Gint,Nint,detF0); 
               
                 Feh_all = zeros(ElementDofs, ElementDofs);

                 H = diag(h * ones(1, ElementDofs));

                 for jj = 1:ElementDofs
                     ukh = uk - H(:,jj);
                     qkh = qk - H(:,jj);
                     Feh_all(:,jj) = Fe_fun(Fe0,ukh,qkh,qk0,phik,Phik,Fibers,Dvec,ElementDofs,PosDofs,Gint,Nint,detF0);
                 end

                 K_loc = (Fe - Feh_all) ./ diag(H)';

           case "AceGen"
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
