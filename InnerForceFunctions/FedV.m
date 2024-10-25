function Fe_dV = FedV(MaterialName,Fe0,u,q,q0,phi,Phi,Fibers,Dvec,Element,ElemDofs,PosDofs,xi,eta,zeta)
ElementName = string(Element);
F_name = 'F_' + ElementName; 
dEde_name = 'dEde_' + ElementName; 
% Exctration necessary parameters 
H = Dvec(end-2); % element' hight   
W = Dvec(end-1); % element' width 
L = Dvec(end);   % element' length
% Adjust the fiber directions 
if Fibers       
   a0 = Dvec(end-6:end-4)';
   a0_name = 'a0_fib_' + ElementName;
   a0_axis = feval(a0_name,a0,q0(PosDofs),phi,Phi,L,H,W,xi,eta,zeta);
   Dvec(end-6:end-4) = a0_axis';
end
const = Dvec(1:end-3);
% Tensor calculations
F = feval(F_name,q,u,q0(PosDofs),phi,L,H,W,xi,eta,zeta);        % Deformation gradient
dEEde = feval(dEde_name,q,u,q0(PosDofs),phi,L,H,W,xi,eta,zeta); % 
SS = feval(MaterialName,F,const);                               %2nd Piola Kirchhoff stress tensor
% Inner force calculations
Fe_dV = Fe0;
for kk=1:ElemDofs  
    for ii=1:3
        for jj=1:3           
            Fe_dV(kk)=Fe_dV(kk)+SS(ii,jj)*dEEde(ii,jj,kk);           
        end
    end
end