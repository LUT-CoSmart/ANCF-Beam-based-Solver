function Fe_dV = FedV(Material,Fe0,ee,Dvec,Element,ElemDofs,xi,eta,zeta)
Fe_dV = Fe0;
const = Dvec(1:end-3);
H = Dvec(end-2); % element hight   
W = Dvec(end-1); % element width 
L = Dvec(end);   % element length

if Element==3243
    F=FF_3243(ee,L,H,W,xi,eta,zeta);
    dEEde=dEde_3243(ee,L,H,W,xi,eta,zeta);
elseif Element==3333
    F=FF_3333(ee,L,H,W,xi,eta,zeta);
    dEEde=dEde_3333(ee,L,H,W,xi,eta,zeta);
elseif Element==3353
    F=FF_3353(ee,L,H,W,xi,eta,zeta);
    dEEde=dEde_3353(ee,L,H,W,xi,eta,zeta);    
elseif Element==3363
    F=FF_3363(ee,L,H,W,xi,eta,zeta);
    dEEde=dEde_3363(ee,L,H,W,xi,eta,zeta);
elseif Element==34103
    F=FF_34X3(ee,L,H,W,xi,eta,zeta);
    dEEde=dEde_34X3(ee,L,H,W,xi,eta,zeta);
end    
%2nd Piola Kirchhoff stress tensor
% Materials. There are Neo-Hookean (1), 2 - and 5 - contants Mooney-Rivlin (2, 5), GOH (0), K.-S. (3).
if Material==1
    SS=Neo(F,const);    
elseif Material==2    
    SS=Mooney2(F,const);
elseif Material==3
    SS=KS(F,const);    
elseif Material==5    
    SS=Mooney5(F,const);
elseif Material==0    
    SS=GOH(F,const);    
end

for kk=1:ElemDofs  
    for ii=1:3
        for jj=1:3           
            Fe_dV(kk)=Fe_dV(kk)+SS(ii,jj)*dEEde(ii,jj,kk);           
        end
   end
end

