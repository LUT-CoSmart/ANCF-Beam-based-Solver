function SS=Neo(F,const)
    mu=const(1);
    d =const(2);
    

    
    II=eye(3);
    J=det(F);
    J_inv23 = J^(-2/3);
    C = F'*F;   % Cauchy-Green deformation tensor
    C_dash = J_inv23*C;
    Cinv = C^(-1);
    C_dash_inv = C_dash^(-1);
    
    
    W_C_dash = mu/2*II;
    
    DEV1 = 0;
    for i = 1:3
        for j = 1:3
            DEV1 = DEV1 + W_C_dash(i,j)*C_dash(j,i);
        end    
    end
    DEV2 = W_C_dash - (1/3)*DEV1*C_dash_inv;
    SS =2/d*(J-1)*J*Cinv+2*J_inv23*DEV2; % The vector of elastic forces without volume integration