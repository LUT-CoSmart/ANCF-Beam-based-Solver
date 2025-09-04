function [Fcont_loc, Ftarg_loc, DOFs_cont, DOFs_targ, Xi_cont, Xi_targ, gap] = Penalty(penalty, ContactBody,TargetBody, Xi)

        gap = abs(Xi(5)); 

        Element_cont = Xi(12);                             % element of slave body 
        DOFs_cont =  ContactBody.xloc(Element_cont,:);     % associated DOFs

        Element_targ = Xi(4);  % element  
        DOFs_targ =  TargetBody.xloc(Element_targ,:);     % associated DOFs    

        Normal = Xi(6:8)'; % NB: it is an original normal (inwards to the respected body)
        % that is why here for Normal_targ and Normal_cont signs are in reverse
                
        Normal_targ = -Normal;
        Normal_cont =  Normal;   

        Fcont_loc =  penalty * gap * Normal_cont;                                                                              
        Ftarg_loc =  penalty * gap * Normal_targ;
    
        Xi_cont = Xi(9:11);
        Xi_targ = Xi(1:4);