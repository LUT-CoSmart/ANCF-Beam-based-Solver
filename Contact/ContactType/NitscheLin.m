function [Fcont_loc, Ftarg_loc, DOFs_cont, DOFs_targ, Xi_cont, Xi_targ, gap] = NitscheLin(penalty, ContactBody,TargetBody, Xi)

        gap = abs(Xi(5)); 
        
        Normal = Xi(6:8)'; % NB: it is an original normal (inwards to the respected body)
        % that is why here for Normal_targ and Normal_cont signs are in reverse

        Normal_targ = -Normal;
        Normal_cont =  Normal;   

        %  Data of slave (contact) body points under the surface of master body 
        xi_cont = Xi(9);
        eta_cont = Xi(10);
        zeta_cont = Xi(11);  
        Xi_cont = Xi(9:11);

        Element_cont = Xi(12);                             % element of slave body 
        DOFs_cont =  ContactBody.xloc(Element_cont,:);     % associated DOFs

        % Data of master (target) body points projected from slave ones
        xi_targ = Xi(1);
        eta_targ = Xi(2);
        zeta_targ = Xi(3);
        Xi_targ = Xi(1:4);

        Element_targ = Xi(4);  % element  
        DOFs_targ =  TargetBody.xloc(Element_targ,:);     % associated DOFs   

        q_targ = TargetBody.q(DOFs_targ);
        u_targ = TargetBody.u(DOFs_targ);
        q0_targ = TargetBody.q0(DOFs_targ);
        phi_targ=TargetBody.phim(Element_targ,:)';
        q0PosDofs_targ = q0_targ(TargetBody.PosDofs);
        F_targ = TargetBody.F(q_targ,u_targ,q0PosDofs_targ,phi_targ,xi_targ,eta_targ,zeta_targ);        % Deformation gradient
        Sigma_targ_n = TargetBody.Sigma_n(F_targ, Normal_targ); 

        q_cont = ContactBody.q(DOFs_cont);
        u_cont = ContactBody.u(DOFs_cont);
        q0_cont = ContactBody.q0(DOFs_cont);
        phi_cont=ContactBody.phim(Element_cont,:)';
        q0PosDofs_cont = q0_cont(ContactBody.PosDofs);
        F_cont = ContactBody.F(q_cont,u_cont,q0PosDofs_cont,phi_cont,xi_cont,eta_cont,zeta_cont);        % Deformation gradient
        Sigma_cont_n = ContactBody.Sigma_n(F_cont, Normal_cont); 
                                     
        % Normal force difference 
        Sigma_n = Sigma_cont_n - Sigma_targ_n;    
        Gap_power = gap;       
        lambda = Gap_power * norm(Sigma_n);
                   
        d_lambda_targ = norm(Sigma_n);
        d_lambda_cont = norm(Sigma_n); 

        Ftarg_loc = (penalty * gap + lambda + Gap_power * d_lambda_targ) * Normal_targ;
        Fcont_loc = (penalty * gap + lambda + Gap_power * d_lambda_cont) * Normal_cont; 