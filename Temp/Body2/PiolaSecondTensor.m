function SS=GOH(F,const)
c10 = const(1);
k1 = const(2);
k2 = const(3);
kappa = const(4);
a0 = const(5:7);

d = const(8);

II=eye(3);
J=det(F);
J_inv23 = J^(-2/3);
C = F'*F;   % Cauchy-Green deformation tensor
C_dash = J_inv23*C; 
Cinv = C^(-1);
C_dash_inv = C_dash^(-1);


A0 = a0' * a0; % a0 is a row-vector

I1=trace(C_dash);
I4=trace(C_dash*A0');  % double contraction operaiton ":";

I4star = kappa*I1 + (1 - 3*kappa)*I4 - 1;
PSI_1 = c10+k1*kappa*exp(k2*I4star^2)*I4star;
PSI_4 = k1*exp(k2*I4star^2)*(1 - 3*kappa)*I4star;


W_C_dash = PSI_1*II + PSI_4*A0;

DEV =  W_C_dash - (1/3)*trace(W_C_dash*C_dash')*C_dash_inv;
SS =2/d*(J-1)*J*Cinv+2*J_inv23*DEV; % The vector of elastic forces without volume integration