function Recovery = RecoveryBeam3333(Body)
        
         h = 1e-3;
         Shape = Body.Shape;
         Gint = Body.Gint; 
         DIM = Body.DIM;         
         Nint = Body.Nint;
               
         % Transform to recovery coordinates
         xi_min   = min(Gint(:,1));
         xi_max   = max(Gint(:,1));
    
         eta_min  = min(Gint(:,2));
         eta_max  = max(Gint(:,2));
    
         zeta_min = min(Gint(:,3));
         zeta_max = max(Gint(:,3));

         N_rec_int = length(xi_interest);
         xi_rec_int = 2*(xi_interest - xi_min)/(xi_max - xi_min) - 1;
         eta_rec_int = 2*(eta_interest - eta_min)/(eta_max - eta_min) - 1;
         zeta_rec_int = 2*(zeta_interest - zeta_min)/(zeta_max - zeta_min) - 1;

         %case 3333,  basis = basis_3333;  required_derivatives = {'y', 'z'};  
         F0_m = F0_points(F0,Gint,Nint); 

         % Stress tensors in each Gauss points 
         Stress_GP   = Stress_points_fun(S,F,Gint,Nint); % quantity itself         
         Stress_GPy  = Stress_points_fun_dir(S,F,F0_m,Gint,Nint,h,2); % first order  derivative y
         Stress_GPz  = Stress_points_fun_dir(S,F,F0_m,Gint,Nint,h,3); % first order  derivative z        

         y = [0;1;0];
         z = [0;0;1];

         Stress_interest =zeros(3,3,N_rec_int);    
         for i = 1:N_rec_int
             xi = xi_rec_int(i);
             eta = eta_rec_int(i);
             zeta = zeta_rec_int(i);                          
             Nm = Shape(xi,eta,zeta); % shape functioit is fucking bullshit!½n values in recovery point
             Nm_reduced = Nm(1,1:DIM:end); % deruced in sense that going from all dofs  
             for j = 1:Nint  

                 for ii = 1:3
                     for jj = 1:3
                         vector = [Stress_GP(ii,jj,j); Stress_GP(ii,jj,j); ]
                         Stress_interest(ii,jj,i) = Nm_reduced *vector;
                     end
                 end    


                 for k = 1:3
                 
                     Stress_interest(:,k,i) = Stress_interest(:,k,i) + Nm_reduced * vector;   
                 end 
             end
         end    

            

end

function F0_m = F0_points(F0,PointSet,Number)
          F0_m = zeros(3,3,Number);
          for k = 1:Number
              xi   = PointSet(k,1);
              eta  = PointSet(k,2);
              zeta = PointSet(k,3);
              F0_m(:,:,k) =  F0(xi,eta,zeta);                  
         end
end          
