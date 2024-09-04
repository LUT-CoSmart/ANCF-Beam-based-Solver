function [K_pp,Fe] = ANCFAceGen(Material,K_pp,Fe,eek,eek0,Dvec,Element,DofsAtNode,Gint,Nint)
addpath('AceGenForces/Element');
uu = eek - eek0;    % Get the vectors of displacements 
% Reshaping to adjust for AceGen
eek0_3 = reshape(eek0, [3, DofsAtNode])';
uu_3 = reshape(uu, [3, DofsAtNode])';
% choosing the suitable AceGen-generated ANCF function
if Element == 3333
    [K_pp,Fe] = ANCFAce3333(Material,eek0_3,uu_3,Dvec,K_pp,Fe,Gint,Nint);    
else

end    
Fe = - Fe; % taking into account the difference between AceGen and Fe_fun          