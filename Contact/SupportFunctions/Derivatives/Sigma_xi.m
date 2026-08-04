function Sigma = Sigma_xi(Body, q, q0, u, xi, eta, zeta)
        
        % _xi means not the derivative by the vectro xi, but dependency of
        % the second Piola tensor from vector xi, it will be used later to
        % obtain the derivatives of the tensor by xi
        F = Body.F(q, q0,u, xi, eta, zeta);
        
        if Body.Fibers            
            S = Body.S(F, a0_fib(xi, eta, zeta) );  
        else
            S = Body.S(F);
        end

        Sigma = 1/det(F) * F * S * F'; 
         
            

            