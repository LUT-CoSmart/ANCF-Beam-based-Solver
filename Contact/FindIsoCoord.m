function xi_eta_zeta_result = FindIsoCoord(Shape,ShapeXi,ShapeEta,ShapeZeta,qk,Point)
    
    point = Point';
    
    Xi = [0; 0; 0];

    r =Shape(Xi(1),Xi(2),Xi(3)) * qk; 

    while norm(r - point) > 1e-7
        
          r_xi = ShapeXi(Xi(1),Xi(2),Xi(3)) * qk;     
          r_eta = ShapeEta(Xi(1),Xi(2),Xi(3)) * qk;     
          r_zeta = ShapeZeta(Xi(1),Xi(2),Xi(3)) * qk;  
          Jac = [r_xi r_eta r_zeta];   

          Xi = Xi - Jac^-1 * (r - point);  
          r =Shape(Xi(1),Xi(2),Xi(3)) * qk;  
    end    
    
    xi_eta_zeta_result = [Xi(1);Xi(2);Xi(3)];