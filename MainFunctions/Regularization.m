function [u,K] = Regularization(K,f,RegType)
        
        if nargin == 2
            RegType = f;
            compute = false;
        elseif nargin == 3
            K_reg = K;
            f_reg = f;
            compute = true;
        else
            error('Incorrect number of inputs in Regularization function');
        end


        switch RegType
               case "off" 
                    
                    
               case {"penalty", "Tikhonov"}
    
                   lambda = 2*sqrt(eps) * norm(K, 'fro');
                   I = eye(size(K));
    
                   if RegType == "penalty"
                      K_reg = K + lambda * I;

                   elseif RegType == "Tikhonov"
                      K_reg = K' * K + lambda * I;
                      
                      if compute
                         f_reg = K' * f;
                      end
                   
                   end    
    
               otherwise
                    error('****** The regularization type is not recognized ******');    
        end
        
        if compute            
            K = sparse(K_reg);
            f = sparse(f_reg);
            % u = - K_reg \ f_reg ;

            tol = 1e-10;
            maxit = 500;
            [u, ~] = pcg(K, f, tol, maxit);
            u = -u;
        else
            u = [];
        end
       


