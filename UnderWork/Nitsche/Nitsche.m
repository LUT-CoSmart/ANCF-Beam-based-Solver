function [Fcont_loc, Ftarg_loc, DOFs_cont, DOFs_targ, Xi_cont, Xi_targ] = Nitsche(penalty, ContactBody,TargetBody, Data)

        Ftarg_loc = zeros(3,1);
        Fcont_loc = zeros(3,1);
        
        Xi_cont = Data.IsoCoordContactPoint;
        Xi_targ = Data.IsoCoorProjected;

        DOFs_cont =  ContactBody.xloc(Data.ElementContactPoint,:); % associated DOFs   
        DOFs_targ =  TargetBody.xloc(Data.ElementProjected,:);

        gap = Data.Gap; % it is negative value

        Normal_targ =  Data.NormalProjected';
        Normal_cont =  Data.ContactPointNormals';      

        % Penalty
        % Ftarg_loc =  penalty * gap * Normal_targ; % normal outward, but gap is negative
        % Fcont_loc =  penalty * gap * Normal_cont; % normal outward, but gap is negative 

        % Nitsche part  
        % h_targ = max(struct2array(TargetBody.Length));
        % nabla_gap = Data.nabla_Gap;
        % nabla_Normal = Data.nabla_Normals;
        % nabla_u_contact = Normal_targ' * nabla_gap + gap*nabla_Normal; 
        % I = eye(3);
        % adding = (gap/h_targ)*nabla_u_contact;   
                
        q_targ= TargetBody.q(DOFs_targ);
        q0_targ= TargetBody.q0(DOFs_targ);
        u_targ= TargetBody.u(DOFs_targ); 
        Phi_targ = TargetBody.Phim(Data.ElementProjected,:); 
        F_targ = TargetBody.F(q_targ,q0_targ,u_targ,Xi_targ(1),Xi_targ(2),Xi_targ(3));         
        S_targ = TargetBody.S(F_targ,q0_targ,Phi_targ,Xi_targ(1),Xi_targ(2),Xi_targ(3));
   
        Sigma_targ = StressRecovery(TargetBody,q_targ,q0_targ,u_targ,Phi_targ,Xi_targ(1),Xi_targ(2),Xi_targ(3));
        %Sigma_targ = TargetBody.Sigma(F_targ,S_targ);

        
        Sigma_targ_nn = Normal_targ' * Sigma_targ * Normal_targ;


        
        q_cont= ContactBody.q(DOFs_cont);
        q0_cont= ContactBody.q0(DOFs_cont); 
        u_cont= ContactBody.u(DOFs_cont);
        Phi_cont = ContactBody.Phim(Data.ElementContactPoint,:);
        F_cont = ContactBody.F(q_cont,q0_cont,u_cont,Xi_cont(1),Xi_cont(2),Xi_cont(3)); 
        S_cont = ContactBody.S(F_cont,q0_cont,Phi_cont,Xi_cont(1),Xi_cont(2),Xi_cont(3));                        
        Sigma_cont = StressRecovery(ContactBody,q_cont,q0_cont,u_cont,Phi_cont,Xi_cont(1),Xi_cont(2),Xi_cont(3));
        Sigma_cont = ContactBody.Sigma(F_cont,S_cont);
        Sigma_cont_nn = Normal_cont' * Sigma_cont * Normal_cont;    
        
        
        % Nitsche contact-pressure trial value.
        trial = Sigma_nn_avg - gap*penalty; 
        if trial > 0 
            Fcont_loc = penalty * gap * Normal_cont - Sigma_nn_avg* Normal_cont; % - nabla_Sigma_nn_cont*gap;
            Ftarg_loc = penalty * gap * Normal_targ - Sigma_nn_avg* Normal_targ ; % - nabla_Sigma_nn_targ*gap; 
        end    

        

