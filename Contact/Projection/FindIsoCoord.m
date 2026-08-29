function xi_eta_zeta_result = FindIsoCoord(Shape,nabla_by_xi,qk,Point)
    
    point = Point';
    
    Xi = [0 0 0];

    r =Shape(Xi(1),Xi(2),Xi(3)) * qk; 

    while norm(r - point) > 1e-5
        
          F_xi = nabla_by_xi(qk,Xi(1),Xi(2),Xi(3)); % 
          Xi = Xi - F_xi^-1 * (r - point);  
          r =Shape(Xi(1),Xi(2),Xi(3)) * qk;  
    end    
    
    xi_eta_zeta_result = [Xi(1) Xi(2) Xi(3)];