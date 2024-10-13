function Fe_dV = FedV(MaterialName,Fe0,u,q,q0,phi,Dvec,Element,ElemDofs,PosDofs,xi,eta,zeta)
Fe_dV = Fe0;
const = Dvec(1:end-3);
H = Dvec(end-2); % element hight   
W = Dvec(end-1); % element width 
L = Dvec(end);   % element length
if Element==3243
    F=F_3243(q,u,q0(PosDofs),phi,L,H,W,xi,eta,zeta);
    dEEde=dEde_3243(q,u,q0(PosDofs),phi,L,H,W,xi,eta,zeta);
elseif Element==3333
    F=F_3333(q,u,q0(PosDofs),phi,L,H,W,xi,eta,zeta);
    dEEde=dEde_3333(q,u,q0(PosDofs),phi,L,H,W,xi,eta,zeta);    
elseif Element==3353
    F=F_3353(q,u,q0(PosDofs),phi,L,H,W,xi,eta,zeta);
    dEEde=dEde_3353(q,u,q0(PosDofs),phi,L,H,W,xi,eta,zeta);    
elseif Element==3363
    F=F_3363(q,u,q0(PosDofs),phi,L,H,W,xi,eta,zeta);
    dEEde=dEde_3363(q,u,q0(PosDofs),phi,L,H,W,xi,eta,zeta);    
elseif Element==34103
    F=F_34103(q,u,q0(PosDofs),phi,L,H,W,xi,eta,zeta);
    dEEde=dEde_34103(q,u,q0(PosDofs),phi,L,H,W,xi,eta,zeta);    
end  
%2nd Piola Kirchhoff stress tensor
SS = feval(MaterialName,F,const);
% Inner force calculations
for kk=1:ElemDofs  
    for ii=1:3
        for jj=1:3           
            Fe_dV(kk)=Fe_dV(kk)+SS(ii,jj)*dEEde(ii,jj,kk);           
        end
   end
end


