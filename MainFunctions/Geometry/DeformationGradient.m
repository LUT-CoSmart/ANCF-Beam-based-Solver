function F = DeformationGradient(SolutionBase,q,q0,u,xi,eta,zeta,nabla_by_xi)

        F0_ = nabla_by_xi(q0,xi,eta,zeta);
        F0_rev = F0_^(-1); 
        I = eye(3);
        
        if SolutionBase == "Position"
           F = nabla_by_xi(q,xi,eta,zeta) * F0_rev;        
        else
           F = I + nabla_by_xi(u,xi,eta,zeta) * F0_rev;
        end
