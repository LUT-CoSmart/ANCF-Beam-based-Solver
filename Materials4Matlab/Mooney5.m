function SS=Mooney5(F,const)
c10 = const(1);
c01 = const(2);
c11 = const(3);
c20 = const(4);
c02 = const(5);
d = const(6);

II=eye(3);
J=det(F);
J_inv23 = J^(-2/3);
C = F'*F;   % Cauchy-Green deformation tensor
C_dash = J_inv23*C;
Cinv = C^(-1);
C_dash_inv = C_dash^(-1);


I1 = trace(C_dash);
I2 = 1/2*(trace(C_dash)^2-trace(C_dash*C_dash'));
W_C_dash = (c10+c11*(I2-3)+2*c20*(I1-3))*II+(c01+c11*(I1-3)+2*c02*(I2-3))*(I1*II-C_dash);
DEV1 = 0;
for i = 1:3
    for j = 1:3
        DEV1 = DEV1 + W_C_dash(i,j)*C_dash(j,i);
    end    
end
DEV2 = W_C_dash - (1/3)*DEV1*C_dash_inv;
SS =2/d*(J-1)*J*Cinv+2*J_inv23*DEV2;