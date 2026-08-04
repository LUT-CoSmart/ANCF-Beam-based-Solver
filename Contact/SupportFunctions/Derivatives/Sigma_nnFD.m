function nabla_Sigma_nn = Sigma_nnFD(Body,N,nabla_r_xi,q,q0,u,xi,eta,zeta)
        
        h = 1e-5; % it is applied in local isoparametric coordinates and variates alongside with element size  

        h_vec = h*eye(3);        
        vec=[xi,eta,zeta]';
        Sigma_xi0 = Sigma_xi(Body,q,q0,u,vec(1),vec(2),vec(3));           
        nabla_Sigma_xi =zeros(9,3);
        
        for i = 1:3
            vec_minus =  vec - h_vec(:,i);

            % For simplicity, stability and consistancy: making sure that the variation of xi is still within the same element
            if vec_minus(i) <-1 
               sign_FD  = +1; 
               sign_FD0 = -1;
               vec_minus =  vec + h_vec(:,i);
            else
               sign_FD  = -1;  
               sign_FD0 = +1;
            end    

            Sigma_xi_minus = Sigma_xi(Body,q,q0,u,vec_minus(1),vec_minus(2),vec_minus(3));
            nabla_Sigma_xi(:,i)= (sign_FD0 * Sigma_xi0(:) + sign_FD * Sigma_xi_minus(:)) / h; 
        end

        nabla_Sigma = nabla_Sigma_xi /nabla_r_xi^(-1);
        nabla_Sigma_nn = zeros(3,1);
        
        for i = 1:3
            nabla_Sigma_nn(i) = N' * reshape(nabla_Sigma(:,i),3,3)  * N;
        end

