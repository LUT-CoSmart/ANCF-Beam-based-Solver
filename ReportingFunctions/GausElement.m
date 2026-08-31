function Set_rec = GausElement(Body,Set)
      
    Set_rec = Set;
    Gint = Body.Gint;        

    % Transform to recovery coordinates
    xi_min   = min(Gint(:,1));
    xi_max   = max(Gint(:,1));

    eta_min  = min(Gint(:,2));
    eta_max  = max(Gint(:,2));

    zeta_min = min(Gint(:,3));
    zeta_max = max(Gint(:,3));
      
    Set_rec(:,1) = 2*(Set(:,1) - xi_min)/(xi_max - xi_min) - 1;
    Set_rec(:,2) = 2*(Set(:,2)- eta_min)/(eta_max - eta_min) - 1;
    Set_rec(:,3) = 2*(Set(:,3)- zeta_min)/(zeta_max - zeta_min) - 1;
