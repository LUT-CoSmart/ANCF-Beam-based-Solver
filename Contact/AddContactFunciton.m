function Body = AddContactFunciton(Body,ContactType)

    L = Body.Length.Ln;
    W = Body.Length.Z;
    H = Body.Length.Y;
    
    Body.F = @(q,u,q0_PosDofs,phi,xi,eta,zeta) F(q,u,q0_PosDofs,phi,L,H,W,xi,eta,zeta);
    Body.nabla_r_xi = @(xi,eta,zeta,q) [Shape_xi_(L,H,W,xi,eta,zeta)*q Shape_eta_(L,H,W,xi,eta,zeta)*q Shape_zeta_(L,H,W,xi,eta,zeta)*q];    
    Body.NodeSphere = feval("MaxNode" + Body.ElementType + "Dimension", Body); % space around node for possible contact check;
    
    if contains(ContactType, "Nitsche")
        % Sigma = @(F_) (1/det(F_) )* F_ * PiolaSecondTensor(F_, Body.const) * F_'; %  Cauchy Stresses 

        Sigma = @(F_) PiolaSecondTensor(F_, Body.const); 
        Body.Sigma_nn = @(F_, N) N'* Sigma(F_) *N;    
        
    end

    if contains(ContactType, "NitscheRigid") || contains(ContactType, "NitscheFull")
        Body.nabla_F = @(q,u,q0_PosDofs,phi,xi,eta,zeta) nabla_F(q,u,q0_PosDofs,phi,L,H,W,xi,eta,zeta);    
        
        % Idea is that we don't know S in advance, so we need to handle it here 
        syms F_syms [3 3] real;
        Sigma_sym = Sigma(F_syms);
        Sigma_F_sym = jacobian(Sigma_sym(:), F_syms(:));

        Sigma_F_symb = matlabFunction(Sigma_F_sym, 'Vars', {F_syms}); % fastes and best option       
        % Sigma_F_auto = @(F_) SigmaJacobian(F_,Sigma);
        % Sigma_F_fd = @(F_) SigmaJacobianFD(F_,Sigma);

        F_xi = @(q,u,q0_PosDofs,phi,xi,eta,zeta) local_nabla_F(q,u,q0_PosDofs,phi,L,H,W,xi,eta,zeta);
        Body.Sigma_F = Sigma_F_symb;
        Body.F_xi = F_xi;
        Body.Sigma_n = @(F_, N) Sigma(F_) * N;
        Body.Sigma = @(F_) Sigma(F_);
        Body.Sigma_xi = @(q,u,q0_PosDofs,phi,xi,eta,zeta) Sigma( F(q,u,q0_PosDofs,phi,L,H,W,xi,eta,zeta) );
    end
