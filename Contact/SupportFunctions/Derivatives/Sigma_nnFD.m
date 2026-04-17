function nabla_Sigma_nn = Sigma_nnFD(N,Sigma_xi,nabla_r_xi,q,u,q0PosDofs,phi,xi,eta,zeta)

        h = 1e-6;
        h_vec = h*eye(3);        
        vec=[xi,eta,zeta]';
        Sigma_xi0 = Sigma_xi(q,u,q0PosDofs,phi,vec(1),vec(2),vec(3));
        nabla_Sigma_xi =zeros(9,3);
        
        for i = 1:3
            %vec_plus = vec + h_vec(:,i);
            vec_minus =  vec - h_vec(:,i);
            %Sigma_xi_plus = Sigma_xi(q,u,q0PosDofs,phi,vec_plus(1),vec_plus(2),vec_plus(3));
            Sigma_xi_minus = Sigma_xi(q,u,q0PosDofs,phi,vec_minus(1),vec_minus(2),vec_minus(3));
            %nabla_Sigma_xi(:,i)= (Sigma_xi_plus(:)-Sigma_xi_minus(:)) / (2*h); 
 
            nabla_Sigma_xi(:,i)= (Sigma_xi0(:)-Sigma_xi_minus(:)) / (h); 
        end

        nabla_Sigma = nabla_Sigma_xi /nabla_r_xi^(-1);
        nabla_Sigma_nn = zeros(3,1);
        for i = 1:3
            nabla_Sigma_nn(i) =  N' * reshape(nabla_Sigma(:,i),3,3)  * N;
        end

