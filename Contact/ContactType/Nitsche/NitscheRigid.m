function [Fcont_loc, Ftarg_loc, DOFs_cont, DOFs_targ, Xi_cont, Xi_targ, gap] = NitscheRigid(penalty, ContactBody,TargetBody, Xi)

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
        Xi_targ = Xi(1:3);

        Element_targ = Xi(4);  % element  
        DOFs_targ =  TargetBody.xloc(Element_targ,:);     % associated DOFs   

        q_targ = TargetBody.q(DOFs_targ);
        u_targ = TargetBody.u(DOFs_targ);
        q0_targ = TargetBody.q0(DOFs_targ);
        phi_targ=TargetBody.phim(Element_targ,:)';
        q0PosDofs_targ = q0_targ(TargetBody.PosDofs);
        F_targ = TargetBody.F(q_targ,u_targ,q0PosDofs_targ,phi_targ,xi_targ,eta_targ,zeta_targ);        % Deformation gradient
        Sigma_targ_nn = TargetBody.Sigma_nn(F_targ, Normal_targ); 


        q_cont = ContactBody.q(DOFs_cont);
        u_cont = ContactBody.u(DOFs_cont);
        q0_cont = ContactBody.q0(DOFs_cont);
        phi_cont=ContactBody.phim(Element_cont,:)';
        q0PosDofs_cont = q0_cont(ContactBody.PosDofs);
        F_cont = ContactBody.F(q_cont,u_cont,q0PosDofs_cont,phi_cont,xi_cont,eta_cont,zeta_cont);        % Deformation gradient
        Sigma_cont_nn = ContactBody.Sigma_nn(F_cont, Normal_cont); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nabla_r_xi_targ = TargetBody.nabla_r_xi(xi_targ,eta_targ,zeta_targ,q_targ);
        nabla_r_xi_cont = ContactBody.nabla_r_xi(xi_cont,eta_cont,zeta_cont,q_cont);

        Sigma_xi_targ = TargetBody.Sigma_xi;
        Sigma_xi_cont = ContactBody.Sigma_xi;
        nabla_Sigma_nn_targ = Sigma_nnFD(Normal_targ,Sigma_xi_targ,nabla_r_xi_targ,q_targ,u_targ,q0PosDofs_targ,phi_targ,xi_targ,eta_targ,zeta_targ);
        nabla_Sigma_nn_cont = Sigma_nnFD(Normal_cont,Sigma_xi_cont,nabla_r_xi_cont,q_cont,u_cont,q0PosDofs_cont,phi_cont,xi_cont,eta_cont,zeta_cont);
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        difference = Sigma_cont_nn - Sigma_targ_nn; % Normal force difference

        lambda = gap * abs(difference); 
        lambda_2_targ =-sign(difference)* gap^2 * nabla_Sigma_nn_targ;
        lambda_2_cont = sign(difference)* gap^2 * nabla_Sigma_nn_cont;

        Ftarg_loc = (penalty * gap + 2*lambda) * Normal_targ + lambda_2_targ;
        Fcont_loc = (penalty * gap + 2*lambda) * Normal_cont + lambda_2_cont; 