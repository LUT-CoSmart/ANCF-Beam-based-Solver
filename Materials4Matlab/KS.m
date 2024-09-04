function SS=KS(F,const)

E=const(1);
nu=const(2);

l=E*nu/((1+nu)*(1-2*nu));
G=E/(2*(1+nu));
%% Saint - Venan material
II=eye(3);
EE=1/2*(F'*F-II);
SS=l*II*trace(EE)+2*G*EE;

