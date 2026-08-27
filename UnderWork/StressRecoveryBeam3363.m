function Stress_interest = StressRecoveryBeam3363(Body,S,F,F0,xi_interest,eta_interest,zeta_interest)
        
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

         %case 3363,  basis = basis_3363;  required_derivatives = {'y', 'z', 'yz', 'yy', 'zz'};   
         F0_m = F0_points(F0,Gint,Nint); 

         % Stress tensors in each Gauss points 
         Stress_GP   = Stress_points_fun(S,F,Gint,Nint); % quantity itself         
         Stress_GPy  = Stress_points_fun_dir(S,F,F0_m,Gint,Nint,h,2); % first order  derivative y
         Stress_GPz  = Stress_points_fun_dir(S,F,F0_m,Gint,Nint,h,3); % first order  derivative z 
         Stress_GPyy = Stress_points_fun_dir2(S,F,F0_m,Gint,Nint,2,2,h);

         Stress_GPzz = Stress_points_fun_dir2(S,F,F0_m,Gint,Nint,3,3,h);
         Stress_GPyz = Stress_points_fun_dir2(S,F,F0_m,Gint,Nint,2,3,h);


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

function Stress_points = Stress_points_fun(S,F,PointSet,Number)                   
         Stress_points = zeros(3,3,Number); % quantity itself
         for k = 1:Number
             xi   = PointSet(k,1);
             eta  = PointSet(k,2);
             zeta = PointSet(k,3);
             F_ = F(xi,eta,zeta);
             Stress_points(:,:,k) =  S(F_,xi,eta,zeta);                  
         end
end

function Stress_points_deriv = Stress_points_fun_dir(S,F,F0_m,PointSet,Number,h,dir)

    % dir = global derivative direction:
    Stress_points_deriv = zeros(3,3,Number);  
    Stress_points = Stress_points_fun(S,F,PointSet,Number);
    
    Stress_points_deriv_local = zeros(3,3,3,Number);% Local derivatives: xi, eta, zeta

    for dir_local = 1:3

        PointSet_h = PointSet;
        PointSet_h(:,dir_local) = PointSet_h(:,dir_local) - h;
        Stress_points_h = Stress_points_fun(S,F,PointSet_h,Number);
        Stress_points_deriv_local(:,:,dir_local,:) = (Stress_points - Stress_points_h)/h;

    end

    for k = 1:Number

        F0_inv = F0_m(:,:,k)^-1;

        Stress_points_deriv(:,:,k) = Stress_points_deriv_local(:,:,1,k) * F0_inv(1,dir) ...
            + Stress_points_deriv_local(:,:,2,k) * F0_inv(2,dir) ...
            + Stress_points_deriv_local(:,:,3,k) * F0_inv(3,dir);

    end

end

function Stress_points_deriv2 = Stress_points_fun_dir2(S,F,F0_m,PointSet,Number,dir1,dir2,h)

    Stress_points_deriv2 = zeros(3,3,Number);

    Stress_local2 = zeros(3,3,3,3,Number);

    Stress_points = Stress_points_fun(S,F,PointSet,Number);


    %Local second 
    for a = 1:3
        for b = 1:3

            if a == b
                PointSet_plus  = PointSet;
                PointSet_minus = PointSet;

                PointSet_plus(:,a)  = PointSet_plus(:,a)  + h;
                PointSet_minus(:,a) = PointSet_minus(:,a) - h;

                Stress_plus = Stress_points_fun(S,F,PointSet_plus,Number);

                Stress_minus = Stress_points_fun(S,F,PointSet_minus,Number);

                Stress_local2(:,:,a,b,:) = (Stress_plus - 2*Stress_points + Stress_minus)/h^2;


            else

                % -----------------------------------------------
                % Mixed derivative:
                %
                % sigma_xieta, sigma_xizeta, sigma_etazeta
                % -----------------------------------------------

                PointSet_pp = PointSet;
                PointSet_pm = PointSet;
                PointSet_mp = PointSet;
                PointSet_mm = PointSet;

                PointSet_pp(:,a) = PointSet_pp(:,a) + h;
                PointSet_pp(:,b) = PointSet_pp(:,b) + h;

                PointSet_pm(:,a) = PointSet_pm(:,a) + h;
                PointSet_pm(:,b) = PointSet_pm(:,b) - h;

                PointSet_mp(:,a) = PointSet_mp(:,a) - h;
                PointSet_mp(:,b) = PointSet_mp(:,b) + h;

                PointSet_mm(:,a) = PointSet_mm(:,a) - h;
                PointSet_mm(:,b) = PointSet_mm(:,b) - h;


                Stress_pp = Stress_points_fun(S,F,PointSet_pp,Number);

                Stress_pm = Stress_points_fun(S,F,PointSet_pm,Number);

                Stress_mp = Stress_points_fun(S,F,PointSet_mp,Number);

                Stress_mm = Stress_points_fun(S,F,PointSet_mm,Number);

                Stress_local2(:,:,a,b,:) = (Stress_pp - Stress_pm - Stress_mp + Stress_mm) /(4*h^2);

            end
        end
    end

    for k = 1:Number

        F0_inv = F0_m(:,:,k)^-1;

        for a = 1:3
            for b = 1:3

                Stress_points_deriv2(:,:,k) = Stress_points_deriv2(:,:,k) + ...
                Stress_local2(:,:,a,b,k) * F0_inv(a,dir1) * F0_inv(b,dir2);

            end
        end

    end

end         


        