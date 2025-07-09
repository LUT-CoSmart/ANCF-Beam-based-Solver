function status = printStatus(deltaf, u_bc, Re, i, ii, imax, steps, titertot, Gap)
    
     persistent previousDisp previousDisp_2

     if i == 1
        previousDisp = 0; 
        previousDisp_2 = 0;
     end
     
     if  nargin < 9 % Gap is optional
         Gap= NaN;
     end   
        
     if  all(abs(deltaf) < Re) || (norm(u_bc)<Re^2) || (abs(2*norm(u_bc) - previousDisp - previousDisp_2)<Re^2) 
                    
         if ~isnan(Gap)
             fprintf('Convergence: %10.4f, Displacements norm: %10.4f, Total gap: %10.7f\n', norm(abs(deltaf)), norm(u_bc), Gap);            
         else
             fprintf('Convergence: %10.4f, Displacements norm: %10.4f\n', norm(abs(deltaf)), norm(u_bc));            
         end
         fprintf('Solution for %d / %d step  is found on %d iteration, Total CPU-time: %f\n', i, steps, ii, titertot);
         status = true;
     elseif ii==imax 
         fprintf('The solution for %d step is not found. The maximum number of iterations is reached. Total CPU-time: %d\n', i, ii);
         status = false;   
     else     
         if ~isnan(Gap)
            fprintf('Iteration: %d, Convergence: %10.4f, Displacements norm: %10.5f, Total gap: %10.7f\n', ii, norm(abs(deltaf)), norm(u_bc), Gap);
         else
            fprintf('Iteration: %d, Convergence: %10.4f, Displacements norm: %10.5f\n', ii, norm(abs(deltaf)), norm(u_bc));
         end
         status = false;
     end 
    
     % Update persistent variables
     if (ii > 1) && (ii < imax)
        previousDisp_2 = previousDisp;
        previousDisp = norm(u_bc);
     else
        previousDisp_2 = 0;
        previousDisp = 0;
     end
    
     

     