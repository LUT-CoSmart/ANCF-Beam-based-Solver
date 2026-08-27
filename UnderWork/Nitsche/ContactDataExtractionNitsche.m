function [Sigma_nn,nabla_Sigma_nn] = ContactDataExtractionNitsche(Body,Element,DOFs,N,Xi)

         xi = Xi(1);
         eta = Xi(2);
         zeta = Xi(3);

         q= Body.q(DOFs);
         u= Body.u(DOFs);
         q0= Body.q0(DOFs);
         Phi = Body.Phim(Element,:);

         F = Body.F(q,q0,u,xi,eta,zeta);    
         S = Body.S(F,q0,Phi,xi,eta,zeta);   

         Sigma = Body.Sigma;

         Sigma_nn =  N' * Sigma(F, S) * N; 
         
         nabla_r_xi = Body.nabla_r_xi(q,xi,eta,zeta);  

         h = 1e-4; % it is applied in local isoparametric coordinates 
                    % and variates alongside with element size (no need being large)  

         nabla_xi_Sigma_nn  = zeros(1,3);
        
         vec=[xi,eta,zeta]';
        
         for i = 1:3

            delta = -h;
            vec_shift =  vec;

            % For simplicity, stability and consistency: 
            % making sure the variation of xi is within the same element            
            if vec(i) - h < -1
                delta = h;
            end

            vec_shift(i) =vec_shift(i) + delta;
            F_shift = Body.F(q,q0,u,vec_shift(1),vec_shift(2),vec_shift(3));    
            S_shift = Body.S(F_shift,q0,Phi,vec_shift(1),vec_shift(2),vec_shift(3));
            Sigma_nn_shift =  N' * Sigma(F_shift,S_shift) * N;

            %% TODO: be sure it in the same patch

            nabla_xi_Sigma_nn(i) = (Sigma_nn_shift - Sigma_nn) / delta;
         end

         nabla_Sigma_nn = nabla_xi_Sigma_nn * nabla_r_xi^(-1);
         nabla_Sigma_nn = nabla_Sigma_nn';
 