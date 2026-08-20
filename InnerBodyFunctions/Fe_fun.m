function Fe = Fe_fun(q,u,F,dF_dq,dEdq,S,ElemDofs,Gint,Nint) 

Fe = zeros(ElemDofs,1);

for ii=1:Nint    % integration all over the element's volume

    xi = Gint(ii,1);
    eta = Gint(ii,2);
    zeta = Gint(ii,3);
    w = Gint(ii,4);
    
    F_ = F(q,u,xi,eta,zeta); % Deformation gradient
    SS = S(F_,xi,eta,zeta);

    % vectorized way
    dF_dq_ = dF_dq(xi,eta,zeta);
    dEdq_ = dEdq(dF_dq_,F_);
    Fe = Fe + reshape(dEdq_,9,ElemDofs).' * SS(:) * w;
    

    % Svec  = SS(:)';
    % dF_dq_ = dF_dq(xi,eta,zeta);
    % 
    % % Inner force calculations
    % for kk=1:ElemDofs 
    %     dF_dq_vec = dF_dq_(:,:,kk);
    %     dEdq_  = dEdq(dF_dq_vec,F_);
    %     Fe(kk) = Fe(kk) + Svec  * dEdq_(:) * w;
    % end    

end

