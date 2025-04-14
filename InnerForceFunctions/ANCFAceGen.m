function [K_loc,Fe] = ANCFAceGen(MaterialName,K_loc,Fe,uk,qk0,qk0f,Dvec,Element,DofsAtNode,Gint,Nint,fac)
% Reshaping to adjust for AceGen
qk0f_3 = reshape(qk0f, [3, DofsAtNode])';
qk0_3 = reshape(qk0, [3, DofsAtNode])';
uk_3 = reshape(uk, [3, DofsAtNode])';
% Construct the function name
function_name = 'ANCF'+string(Element)+MaterialName;
[~,~,~,~,K_loc,Fe,~,~] = feval(function_name,qk0f_3,qk0_3,uk_3,Dvec,K_loc,Fe,Gint',Nint);
K_loc = fac*K_loc;
Fe = - fac*Fe; % taking into account the difference between AceGen and Fe_fun          