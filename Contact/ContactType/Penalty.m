function [Fcont_loc, Ftarg_loc] = Penalty(penalty, ContactBody,TargetBody, Data)

        gap = abs(Data.Gap); 

        Normal = Data.Normal; 
                
        Normal_targ =  Normal;
        Normal_cont = -Normal;   
       
        Fcont_loc =  penalty * gap * Normal_cont;                                                                              
        Ftarg_loc =  penalty * gap * Normal_targ;
  