function [Recovery_S1,Recovery_S2] = RecoveryMain(Body)
 
    Shape_        = Body.Shape;
    xloc          = Body.xloc;
    TotalDofs     = Body.TotalDofs;
    ElementNumber = Body.ElementNumber;
    q             = Body.q;
    q0            = Body.q0;
    q0f           = Body.q0f;
    u             = Body.u;     
    Gint          = Body.Gint;
    Nint          = Body.Nint;
    Sigma         = Body.Sigma;
    ElementDofs   = Body.ElementDofs;
    ElemDim       = 0.5 * Body.Length.Ln * Body.detF0; 
    
    M = zeros(TotalDofs, TotalDofs);
    b = zeros(TotalDofs, 1);

    b_S1 = zeros(TotalDofs, 1);
    b_S2 = zeros(TotalDofs, 1);
    for ii = 1:ElementNumber

        edofs = xloc(ii, :);
        qk    = q(edofs);
        qk0   = q0(edofs);
        qk0f  = q0f(edofs);
        uk    = u(edofs);

        Mlocal = zeros(ElementDofs);
        blocal_q = zeros(ElementDofs, 1);

        blocal_S1 = zeros(ElementDofs, 1);
        blocal_S2 = zeros(ElementDofs, 1);
        
        for j = 1:Nint

            xi = Gint(j, 1);
            eta = Gint(j, 2);
            zeta = Gint(j, 3);
            w = Gint(j, 4);

            N = Shape_(xi, eta, zeta);            
            Mlocal = Mlocal + N' * N * w * ElemDim;

            % for q 
            r = N * qk;
            blocal_q = blocal_q + N' * r * w * ElemDim;
            
            % for q_f
            F_ = Body.F(qk,qk0,uk,xi,eta,zeta);
            S =  Body.S(F_,qk0,qk0f,xi,eta,zeta); 
            S_m = Sigma(F_,S);
            Sigma_vec = [S_m(1,1) S_m(2,2) S_m(3,3) S_m(1,2) S_m(2,3) S_m(1,3)]'; 

            blocal_S1 = blocal_S1 + N' * Sigma_vec(1:3) * w * ElemDim;
            blocal_S2 = blocal_S2 + N' * Sigma_vec(4:6) * w * ElemDim;
        end

        M(edofs, edofs) = M(edofs, edofs) + Mlocal;
        b(edofs)        = b(edofs)        + blocal_q;
        b_S1(edofs)     = b_S1(edofs)     + blocal_S1;
        b_S2(edofs)     = b_S2(edofs)     + blocal_S2;
    end

    Recovery_q  = M \ b;
    fprintf('q relative recovery error: %.3e\n', norm(q - Recovery_q) / max(norm(q), eps));

    Recovery_S1  = M \ b_S1;
    Recovery_S2  = M \ b_S2;

    

 