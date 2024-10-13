function [K,ff,Fext]=Kt(q,q0,u,phim,F,xloc,nx,nn,nl,MaterialName,Dvec,detF0,Gint,Nint,sol_acegen,Element,ElemDofs,DofsAtNode,PosDofs,h)
addpath('InnerForceFunctions');
Fext=Fext_ANCF(DofsAtNode,nn,nx,F);
K=zeros(nx,nx);
Fint=zeros(nx,1); 
for ii = 1:nl
    uk=u(xloc(ii,:));        
    qk=q(xloc(ii,:));    
    qk0=q0(xloc(ii,:)); 
    phik=phim(ii,:)';    
    K_loc = zeros(ElemDofs,ElemDofs);  % initialization of local stiffness matrix
    Fe   = zeros(ElemDofs,1);          % initialization of local inner force vector
    if sol_acegen == false        
        [K_loc,Fe] = ANCFMatlab(MaterialName,K_loc,Fe,uk,qk,qk0,phik,Dvec,Element,ElemDofs,PosDofs,Gint,Nint,detF0,h); 
    else
        [K_loc,Fe] = ANCFAceGen(MaterialName,K_loc,Fe,uk,qk0,Dvec,Element,DofsAtNode,Gint,Nint);
    end       
    for jj = 1:ElemDofs
    	ind01 = xloc(ii,jj); %Index 01
        for kk = 1:ElemDofs
            ind02 = xloc(ii,kk); % Index 02
            K(ind01,ind02) = K(ind01,ind02)+K_loc(jj,kk);
        end
        Fint(ind01) = Fint(ind01)+Fe(jj); %%% sign infron of F
     end %End of Assembly pp
end    
ff = Fint - Fext;