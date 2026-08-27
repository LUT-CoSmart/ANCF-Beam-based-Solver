function [Fcont_loc, Ftarg_loc, DOFs_cont, DOFs_targ, Xi_cont, Xi_targ] = Penalty(penalty, ContactBody,TargetBody, Data)
        
        Xi_cont = Data.IsoCoordContactPoint;
        Xi_targ = Data.IsoCoorProjected;

        DOFs_cont =  ContactBody.xloc(Data.ElementContactPoint,:); % associated DOFs   
        DOFs_targ =  TargetBody.xloc(Data.ElementProjected,:);

        gap = Data.Gap; 

        Normal_targ =  Data.NormalProjected';
        Normal_cont =  Data.ContactPointNormals';   

        Ftarg_loc =  penalty * gap * Normal_targ; % normal outward, but gap is negative
        Fcont_loc =  penalty * gap * Normal_cont; % normal outward, but gap is negative                                                                             
        
