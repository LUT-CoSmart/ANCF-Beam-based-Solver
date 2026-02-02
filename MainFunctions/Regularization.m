function u = Regularization(K,f,RegType)
        
      if RegType ~= "off"                
                x0 = zeros(length(f_reg),1);
                [u, ~] = pcg(K, f, sqrt(eps), 500, [], [], x0);
                u = -u;                
      else
            num = cond(K);
            if num > 1e12
               D = diag(1./sqrt(sum(K.^2,2)));
               u =- (D * K)\(D * f);
            else 
               u = -K\f;
            end
      end  
     



