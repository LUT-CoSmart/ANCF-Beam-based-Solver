function b_f = RecoveryMain_v2(Body)
 
    xloc          = Body.xloc;
    TotalDofs     = Body.TotalDofs;
    ElementNumber = Body.ElementNumber;
    q             = Body.q;
    q0            = Body.q0;
    q0f           = Body.q0f;
    u             = Body.u;     
    Gint          = Body.Gint;
    Nint          = Body.Nint;
    ElementDofs   = Body.ElementDofs;
    ElemDim       = 0.5 * Body.Length.Ln * Body.detF0; 
    dEdq          = Body.dEdq;

    b_f = zeros(TotalDofs, 1);
    b_f_or = zeros(TotalDofs, 1);
    b_dir = zeros(TotalDofs, 1);
    for ii = 1:ElementNumber

        edofs = xloc(ii, :);
        qk    = q(edofs);
        qk0   = q0(edofs);
        qk0f  = q0f(edofs);
        uk    = u(edofs);
    
        dF_dq = @(xi,eta,zeta) Body.dF_dq(qk0,xi,eta,zeta);
        S = @(F_,xi,eta,zeta) Body.S(F_,qk0,qk0f,xi,eta,zeta);  
        F = @(q,u,xi,eta,zeta) Body.F(q,qk0,u,xi,eta,zeta);

        fint_local = zeros(ElementDofs,1);
        fint_local_or = zeros(ElementDofs,1);
        
        fint_local_dir = Fe_fun(qk,uk,F,dF_dq,dEdq,S,ElementDofs,Gint,Nint) * ElemDim;
        

        for j = 1:Nint
        
            xi   = Gint(j,1);
            eta  = Gint(j,2);
            zeta = Gint(j,3);
            w    = Gint(j,4);
        
            F_ = F(qk,uk,xi,eta,zeta);            
            S_ = S(F_,xi,eta,zeta); % Second Piola–Kirchhoff stress
             
            % First Piola–Kirchhoff stress, conjugate to F
            P_ = F_ * S_;
            dF_dq_ = dF_dq(xi,eta,zeta);
            dEdq_ = dEdq(dF_dq_,F_);

            D_F = reshape(dF_dq_,9,ElementDofs);
            D_E = reshape(dEdq_,9,ElementDofs);

            fint_local    = fint_local    + D_F.' * P_(:) * w * ElemDim;
            fint_local_or = fint_local_or + D_E.' * S_(:) * w * ElemDim; 
            
        end
        b_f(edofs) = b_f(edofs) + fint_local;
        b_f_or(edofs) = b_f_or(edofs) + fint_local_or;
        b_dir(edofs) = b_dir(edofs) + fint_local_dir;
    end

den = max(norm(Body.Fe),eps);

err = norm(b_f - Body.Fe) / den;
err_or  = norm(b_f_or - Body.Fe) / den;
err_dir  = norm(b_dir - Body.Fe) / den;

fprintf('||b_f - Fe||/||Fe|| = %.3e\n',err);
fprintf('||b_f_or - Fe||/||Fe|| = %.3e\n',err_or);
fprintf('||b_dir - Fe||/||Fe|| = %.3e\n', err_dir);
