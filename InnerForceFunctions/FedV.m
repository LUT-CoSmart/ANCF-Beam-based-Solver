function Fe_dV = FedV(Fe0,u,q,q0,phi,Phi,Fibers,Dvec,ElemDofs,PosDofs,xi,eta,zeta)

    dEde_ = zeros(3,3,ElemDofs);
    SS = zeros(3,3);
    % Exctration necessary parameters 
    H = Dvec(end-2); % element' hight   
    W = Dvec(end-1); % element' width 
    L = Dvec(end);   % element' length
    % Adjust the fiber directions 
    if Fibers       
       a0 = Dvec(end-6:end-4)';
       a0_axis = a0_fib(a0,q0(PosDofs),phi,Phi,L,H,W,xi,eta,zeta);
       Dvec(end-6:end-4) = a0_axis';
    end
    const = Dvec(1:end-3);
    % Tensor calculations
    F_ = F(q,u,q0(PosDofs),phi,L,H,W,xi,eta,zeta);        % Deformation gradient
    dEde_ = dEde(q,u,q0(PosDofs),phi,L,H,W,xi,eta,zeta);  %        
    SS = PiolaSecondTensor(F_, const);                        

    % Inner force calculations
    Fe_dV = Fe0;
    for kk=1:ElemDofs  
        for ii=1:3
            for jj=1:3           
                Fe_dV(kk)=Fe_dV(kk)+SS(ii,jj)*dEde_(ii,jj,kk);           
            end
        end
    end