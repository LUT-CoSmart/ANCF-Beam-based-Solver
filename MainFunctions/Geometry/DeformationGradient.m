function F = DeformationGradient(SolutionBase,q,q0,u,L,H,W,xi,eta,zeta)

        F0_ = F0(q0,L,H,W,xi,eta,zeta);
        F0_rev = F0_^(-1); 
        I = eye(3);
        
        if SolutionBase == "Position"
           F = F_xi(q,L,H,W,xi,eta,zeta) * F0_rev;        
        else
           F = I + F_xi(u,L,H,W,xi,eta,zeta) * F0_rev;
        end
