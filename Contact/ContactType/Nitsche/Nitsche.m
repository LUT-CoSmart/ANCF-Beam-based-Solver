function [Fcont_loc, Ftarg_loc, DOFs_cont, DOFs_targ, Xi_cont, Xi_targ, gap] = Nitsche(penalty, ContactBody,TargetBody, Xi)

        gap = abs(Xi(5)); 
        
        Xi_cont = Xi(9:11);
        Xi_targ = Xi(1:3);

        Element_cont = Xi(12);                             % element of slave body 
        DOFs_cont =  ContactBody.xloc(Element_cont,:);     % associated DOFs

        Element_targ = Xi(4);  % element  
        DOFs_targ =  TargetBody.xloc(Element_targ,:);     % associated DOFs    

        Normal = Xi(6:8)'; % NB: it is an original normal (inwards to the respected body)
                           % that is why here for Normal_targ and Normal_cont signs are in reverse
                
        Normal_targ = -Normal;
        Normal_cont =  Normal;   
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [Sigma_nn_cont,nabla_Sigma_nn_cont] = ContactDataExtractionNitsche(ContactBody,Element_cont,DOFs_cont,Normal_cont,Xi_cont);
        [Sigma_nn_targ,nabla_Sigma_nn_targ] = ContactDataExtractionNitsche(TargetBody,Element_targ,DOFs_targ,Normal_targ,Xi_targ);
       
        Sigma_nn_avg = 0.5 * (Sigma_nn_cont + Sigma_nn_targ);

        % Nitsche contact-pressure trial value.
        trial = gap - Sigma_nn_avg/penalty - (nabla_Sigma_nn_cont'*Normal_cont)/(2*penalty)
    
        active = trial > 0;
    
        if active
            [ Sigma_nn_avg nabla_Sigma_nn_cont'*Normal_cont]

           Fcont_loc = penalty * trial * Normal_cont - nabla_Sigma_nn_cont*gap;
           Ftarg_loc = penalty * trial * Normal_cont - nabla_Sigma_nn_targ*gap; 
        else
           Fcont_loc = zeros(3,1);
           Ftarg_loc = zeros(3,1);
        end

        
end