function [K_pp,Fe] = ANCFAce3333(Material,eek0_3,uu_3,Dvec,K_pp,Fe,Gint,Nint)
    addpath('AceGenForces/3333');    
    if Material == 1 
            [~,~,~,K_pp,Fe,~,~] = ANCF3333Neocont(eek0_3,uu_3,Dvec,K_pp,Fe,Gint',Nint);
    elseif Material == 3
            [~,~,~,K_pp,Fe,~,~] = ANCF3333KScont(eek0_3,uu_3,Dvec,K_pp,Fe,Gint',Nint);
    end
    
    