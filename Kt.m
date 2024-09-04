function [K,ff,Fext]=Kt(ee,ee0,F,xloc,nx,nn,nl,Material,Dvec,detF0,Gint,Nint,sol,Element,ElemDofs,DofsAtNode,h)
addpath('ForceFunctions');
Fext=Fext_ANCF(DofsAtNode,nn,nx,F);
K=zeros(nx,nx);
Fint=zeros(nx,1); 
for ii = 1:nl
    eek=ee(xloc(ii,:));    
    eek0=ee0(xloc(ii,:));
    K_pp = zeros(ElemDofs,ElemDofs);  % initialization of local stiffness matrix
    Fe   = zeros(ElemDofs,1);         % initialization of local inner force vector
    if sol == 0
        [K_pp,Fe] = ANCFMatlab(Material,K_pp,Fe,eek,Dvec,Element,ElemDofs,Gint,Nint,detF0,h); 
    else
        [K_pp,Fe] = ANCFAceGen(Material,K_pp,Fe,eek,eek0,Dvec,Element,DofsAtNode,Gint,Nint);
    end       
    for jj = 1:ElemDofs
    	ind01 = xloc(ii,jj); %Index 01
        for kk = 1:ElemDofs
            ind02 = xloc(ii,kk); % Index 02
            K(ind01,ind02) = K(ind01,ind02)+K_pp(jj,kk);
        end
        Fint(ind01) = Fint(ind01)+Fe(jj); %%% sign infron of F
     end %End of Assembly pp
end    
ff = Fint - Fext;