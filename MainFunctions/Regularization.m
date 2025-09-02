function [u,K_reg] = Regularization(K,f,RegType,compute)
        
        if nargin == 3
            compute = true;
        end
        
        if  nargin < 3   
            error('Incorrect number of inputs in Regularization function');
        end


        switch RegType
               case "off" 
                  K_reg = K;
                  f_reg = f;  
                    
               case {"penaltyKf", "penaltyK", "Tikhonov", "Papadopoulos" }
                                         
                   I = eye(size(K,2));

                   if RegType == "penaltyKf"
                       

                      lambda = sqrt(eps)  * norm(K, 'fro') / max(norm(f),1);
                      K_reg = K + lambda * I;
                      f_reg = f;
                    
                   elseif RegType == "penaltyK"
                       
                      lambda = sqrt(eps) * norm(K, 'fro'); 
                      K_reg = K + lambda * I;
                      f_reg = f;

                   elseif RegType == "Tikhonov"
                      lambda = sqrt(eps) * norm(K, 'fro'); 
                      K_reg = K' * K + lambda * I;
                      f_reg = K' * f;
                   elseif RegType == "Papadopoulos"
                      lambda = sqrt(eps) * norm(K, 'fro'); 
                      K_reg = K + lambda * I;
                      f_reg = f + lambda * ones(length(f),1);

                  end    
    

               otherwise
                    error('****** The regularization type is not recognized ******');    
        end
        
        if compute  
             
            % if ~issparse(K_reg), K_reg = sparse(K_reg); end
            % if ~issparse(f_reg), f_reg = sparse(f_reg); end

            if RegType ~= "off"                
    
                tol = sqrt(eps);
                maxit = 500;
    
                x0 = zeros(length(f_reg),1);
                [u, ~] = pcg(K_reg, f_reg, tol, maxit, [], [], x0);
                u = -u;
                
            else

                u = - K_reg\ f_reg;

            end   

        else
            u = [];
        end
       


