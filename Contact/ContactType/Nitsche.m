function [Fcont_loc, Ftarg_loc, DOFs_cont, DOFs_targ, Xi_cont, Xi_targ, gap] = Nitsche(ContactType,penalty, ContactBody,TargetBody, Xi)

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
        
        Element_targ = Xi(4);  % element  
        DOFs_targ =  TargetBody.xloc(Element_targ,:);     % associated DOFs   
        
        % Data of master (target) body points projected from slave ones
        xi_targ = Xi(1);
        eta_targ = Xi(2);
        zeta_targ = Xi(3);
        Xi_targ = Xi(1:3);

        q_targ = TargetBody.q(DOFs_targ);
        u_targ = TargetBody.u(DOFs_targ);
        q0_targ = TargetBody.q0(DOFs_targ);
        F_targ = TargetBody.F(q_targ,q0_targ,u_targ,xi_targ,eta_targ,zeta_targ); 
               
        % taking into account fiber enriched material models 
        if TargetBody.Fibers
            Phi_targ = TargetBody.Phim(Element_targ,:);
            a0_axis_targ =  TargetBody.a0_fib(q0_targ,Phi_targ,xi_targ,eta_targ,zeta_targ);
            S_targ = TargetBody.S(F_targ, a0_axis_targ);  
        else
            S_targ = TargetBody.S(F_targ);  
        end        
        Sigma_targ_nn = TargetBody.Sigma_nn(F_targ, Normal_targ, S_targ); 

        q_cont = ContactBody.q(DOFs_cont);
        u_cont = ContactBody.u(DOFs_cont);
        q0_cont = ContactBody.q0(DOFs_cont);
        F_cont = ContactBody.F(q_cont,q0_cont,u_cont,xi_cont,eta_cont,zeta_cont); 
        
        % taking into account fiber enriched material models 
        if ContactBody.Fibers
            Phi_cont = ContactBody.Phim(Element_cont,:);
            a0_axis_cont =  ContactBody.a0_fib(q0_cont,Phi_cont,xi_cont,eta_cont,zeta_cont);
            S_cont = ContactBody.S(F_cont, a0_axis_cont);  
        else
            S_cont = ContactBody.S(F_cont);  
        end                
        Sigma_cont_nn = ContactBody.Sigma_nn(F_cont, Normal_targ, S_cont); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Normal force difference
        difference = (Sigma_cont_nn - Sigma_targ_nn); % Normal force difference

        % Stabilization
        alpha = min(1, 1.5 * penalty/ (norm(difference)+1e-7) ); % No more than 150% of penalty
                                                                 % Small eps to avoid zero division                 
        % alpha = 1;
        lambda =  alpha *  norm(difference); % added division on penalty for the regularization
        %lambda = 0;

        Fcont_loc = gap * (penalty + 2*lambda) * Normal_cont; 
        Ftarg_loc = gap * (penalty + 2*lambda) * Normal_targ; % + due to different normal orientation

        % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        if ContactType ~= "NitscheLin"

           Multiplier = alpha * sign(difference) * gap^2; % added division on penalty for the regularization
            
           nabla_r_xi_targ = TargetBody.nabla_r_xi(xi_targ,eta_targ,zeta_targ,q_targ);
           nabla_r_xi_cont = ContactBody.nabla_r_xi(xi_cont,eta_cont,zeta_cont,q_cont);
           
           % here Normal doesn't matter, because we use (n ⊗ n)
           % Attention to the sign, such as we use it for pressure         
           nabla_Sigma_nn_targ = -Sigma_nnFD(TargetBody,Normal_targ,nabla_r_xi_targ,q_targ,q0_targ,u_targ,xi_targ,eta_targ,zeta_targ);
           nabla_Sigma_nn_cont = -Sigma_nnFD(ContactBody,Normal_cont,nabla_r_xi_cont,q_cont,q0_cont,u_cont,xi_cont,eta_cont,zeta_cont);

           lambda_2_targ =-Multiplier * nabla_Sigma_nn_targ; % due to different normal orientation
           lambda_2_cont = Multiplier * nabla_Sigma_nn_cont; 

           Ftarg_loc = Ftarg_loc + lambda_2_targ;
           Fcont_loc = Fcont_loc + lambda_2_cont; 
        end    
