function [K_loc,Fe] = ANCFAceGen(MaterialName,K_loc,Fe,uk,qk0,Dvec,Element,DofsAtNode,Gint,Nint)
% Go to the required folder 
path = 'TensorDerivations/AceGen/'+string(Element);
addpath(path);
% Reshaping to adjust for AceGen
qk0_3 = reshape(qk0, [3, DofsAtNode])';
uk_3 = reshape(uk, [3, DofsAtNode])';
% Construct the function name
function_name = 'ANCF'+string(Element)+MaterialName;
[~,~,~,K_loc,Fe,~,~] = feval(function_name,qk0_3,uk_3,Dvec,K_loc,Fe,Gint',Nint);
Fe = - Fe; % taking into account the difference between AceGen and Fe_fun          