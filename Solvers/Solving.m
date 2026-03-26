function u = Solving(K,f)
        
            num = cond(K);
            if num > 1e12
                rowNorms = sqrt(sum(K.^2, 2));
                rowNorms(rowNorms == 0) = 1;
                D = diag(1 ./ rowNorms);
                u =- (D * K)\(D * f);
            else 
               u = -(K\f);
            end

      
   



