function Fe = Fe_fun(q,u,F,dF_dq,dEdq,S,ElemDofs,Gint,Nint) 

Fe = zeros(ElemDofs,1);

for ii=1:Nint    % integration all over the element's volume

    xi = Gint(ii,1);
    eta = Gint(ii,2);
    zeta = Gint(ii,3);
    w = Gint(ii,4);
    
    F_ = F(q,u,xi,eta,zeta); % Deformation gradient
    SS = S(F_,xi,eta,zeta);

    dF_dq_ = dF_dq(xi,eta,zeta);
    dEdq_ = dEdq(dF_dq_,F_);
    Fe = Fe + reshape(dEdq_,9,ElemDofs).' * SS(:) * w;
 
end

