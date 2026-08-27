function Stress = StressRecovery(Body,F,q0,Phi,xi_interest,eta_interest,zeta_interest)
      
    DIM = Body.DIM;  
    ShapeFunctionNumber = Body.ElementDofs/DIM;
    Shape = Body.Shape;
    Gint = Body.Gint;        
    S = @(F_,xi,eta,zeta) Body.S(F_,q0,Phi,xi,eta,zeta);  

    % Transform to recovery coordinates
    xi_min   = min(Gint(:,1));
    xi_max   = max(Gint(:,1));

    eta_min  = min(Gint(:,2));
    eta_max  = max(Gint(:,2));

    zeta_min = min(Gint(:,3));
    zeta_max = max(Gint(:,3));

    N_rec = Body.Nint;
    A_rec = zeros(N_rec,ShapeFunctionNumber); % first matrix
    xi_rec = 2*(Gint(:,1) - xi_min)/(xi_max - xi_min) - 1;
    eta_rec = 2*(Gint(:,2) - eta_min)/(eta_max - eta_min) - 1;
    zeta_rec = 2*(Gint(:,3) - zeta_min)/(zeta_max - zeta_min) - 1;
    
    % info collecting     
    Stress_GP = zeros(3,3,N_rec);
    for k = 1:N_rec
        xi   = Gint(k,1);
        eta  = Gint(k,2);
        zeta = Gint(k,3);
        F_ = F(xi,eta,zeta);    
        S_ = S(F_,xi,eta,zeta);       
        Stress_GP(:,:,k) = Body.Sigma(F_,S_);
        Nm_rec = Shape(xi_rec(k),eta_rec(k),zeta_rec(k));
        A_rec(k,:) = Nm_rec(1,1:DIM:end);
    end
    
    N_rec_int = length(xi_interest);
    A_rec_int = zeros(N_rec_int,ShapeFunctionNumber);
    xi_rec_int = 2*(xi_interest - xi_min)/(xi_max - xi_min) - 1;
    eta_rec_int = 2*(eta_interest - eta_min)/(eta_max - eta_min) - 1;
    zeta_rec_int = 2*(zeta_interest - zeta_min)/(zeta_max - zeta_min) - 1;

    for k = 1:N_rec_int
        Nm_rec = Shape(xi_rec_int(k),eta_rec_int(k),zeta_rec_int(k));
        A_rec_int(k,:) = Nm_rec(1,1:DIM:end);
    end
    
    Stress = zeros(3,3,N_rec_int);
    for i = 1:3
        for j = 1:3
            stress_all_points = Stress_GP(i,j,:);
            stress_all_points = stress_all_points(:);
            coeff = A_rec \ stress_all_points;
            Stress(i,j,:) = A_rec_int * coeff;
        end
    end
