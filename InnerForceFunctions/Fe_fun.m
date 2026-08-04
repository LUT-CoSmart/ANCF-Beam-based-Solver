function Fe = Fe_fun(Nm_xi,q,u,ElemDim,F,F0,S,ElemDofs,Gint,Nint,DeformationType,SolutionBase) 

Fe = zeros(ElemDofs,1);

for ii=1:Nint    % integration all over the element's volume

    xi = Gint(ii,1);
    eta = Gint(ii,2);
    zeta = Gint(ii,3);
    w = Gint(ii,4);
    
    F0_ = F0(xi,eta,zeta);
    F0_rev = F0_^(-1);
    
    F_ = F(q,u,xi,eta,zeta); % Deformation gradient

    Nm_xi_ = Nm_xi(xi,eta,zeta);

    SS = S(F_,xi,eta,zeta);  

    Svec  = SS(:)';

    % Inner force calculations
    for kk=1:ElemDofs 

        dF_dq_vec = Nm_xi_(:,kk);
        dF_dq_ = reshape(dF_dq_vec, 3, 3) * F0_rev;

        if DeformationType == "Small"
           dEdq = 0.5 * ( dF_dq_' + dF_dq_ );  
        else
           dEdq = 0.5 * (dF_dq_' * F_ + F_' * dF_dq_ );
        end    
        
        Fe(kk) = Fe(kk) + Svec  * dEdq(:) * w;

    end    

end
Fe=Fe*ElemDim;
