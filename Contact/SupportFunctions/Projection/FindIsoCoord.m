function xi_eta_zeta_result = FindIsoCoord(Shape,nabla_r_xi,qk,Point)
    
    point = Point';
    
    Xi = [0; 0; 0];

    r =Shape(Xi(1),Xi(2),Xi(3)) * qk; 

    while norm(r - point) > 1e-7
        
          Jac = nabla_r_xi(Xi(1),Xi(2),Xi(3),qk);
          Xi = Xi - Jac^-1 * (r - point);  
          r =Shape(Xi(1),Xi(2),Xi(3)) * qk;  
    end    
    
    xi_eta_zeta_result = [Xi(1);Xi(2);Xi(3)];