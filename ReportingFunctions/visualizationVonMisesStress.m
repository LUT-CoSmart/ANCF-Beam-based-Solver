function visualizationVonMisesStress(Body,Show)
    
    
    if Show 
        figure();
        axis equal
        visualization(Body,Body.q,'none',Show); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        Shape_ = Body.Shape;
        ElementNumber = Body.ElementNumber;
        Gint = Body.Gint_vis;
        Nint = Body.Nint_vis;
        
        points = [];
        for i = 1:ElementNumber  % parfor doesn't work yet, so == for  
                     
                xloc = Body.xloc;
                uk=Body.u(xloc(i,:));  
                qk=Body.q(xloc(i,:)); 
                qk0=Body.q0(xloc(i,:)); 
                Phik=Body.Phim(i,:)';
                F = @(xi,eta,zeta) Body.F(qk,qk0,uk,xi,eta,zeta);
                
                for ii=1:Nint
                    xi = Gint(ii,1);
                    eta = Gint(ii,2);
                    zeta = Gint(ii,3);
                    
                    F_ = F(xi,eta,zeta); % Deformation gradient    
                    S = Body.S(F_,qk0,Phik,xi,eta,zeta);  % 2nd Piola stress tensor
                    Sigma = Body.Sigma(F_,S); % Cauchy stress tensor 
                    Sigma_dev = Sigma - trace(Sigma)/3 * eye(3); % Deviatoric stress tensor
                    Sigma_VM = sqrt(3/2 * sum(Sigma_dev.^2, 'all')); % von Mises stress
                    
                    r = Shape_(xi,eta,zeta)*qk; 

                    points = [points; r' Sigma_VM];
                end

        end
        h = scatter3(points(:,1), points(:,2), points(:,3), 25, points(:,4), 'filled','MarkerEdgeColor', 'k');        
        h.DataTipTemplate.DataTipRows(4).Label = '\sigma';
        h.DataTipTemplate.DataTipRows(4).Format = '%.6e';

        colormap(turbo);
        cb = colorbar;
        cb.Label.String = 'von Mises stress';
        title('Von Mises stress');

    else
        disp('No visualization')
    end   
    