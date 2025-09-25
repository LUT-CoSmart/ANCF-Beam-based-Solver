function status = printStatus(deltaf, u_bc, Re, i, ii, imax, steps, titertot, Gap)
    
     persistent previousForce previousDisp

     if ii == 1
        previousForce = 0; 
        previousDisp = 0; 
     end
     
     if  nargin < 9 % Gap is optional
         Gap= NaN;
     end   
        
      if  all(abs(deltaf) < Re) %    || ... % standard condition of exit
         %  (norm(u_bc)<Re^2)  || abs(norm(u_bc) - previousDisp)<Re^2 || ... % stop when the change is too small
         % ( all(abs(deltaf) - previousForce)<Re && abs(norm(u_bc) - previousDisp)<Re ) % additional condition for exit  

         if ~isnan(Gap)
             fprintf('Convergence: %10.4f, Displacements norm: %10.4f, Total gap: %10.7f\n', norm(abs(deltaf)), norm(u_bc), Gap);            
         else
             fprintf('Convergence: %10.4f, Displacements norm: %10.4f\n', norm(abs(deltaf)), norm(u_bc));            
         end
         fprintf('Solution for %d / %d step  is found on %d iteration, Total CPU-time: %.2f\n', i, steps, ii, titertot);
         status = true;
     elseif ii==imax 
         fprintf('The solution for %d step is not found. The maximum number of iterations is reached. Total CPU-time: %.2f\n', i, titertot);
         status = false;   
     else     
         if ~isnan(Gap)
            fprintf('Iteration: %d, Convergence: %10.4f, Displacements norm: %10.5f, Total gap: %10.7f\n', ii, norm(abs(deltaf)), norm(u_bc), Gap);
         else
            fprintf('Iteration: %d, Convergence: %10.4f, Displacements norm: %10.5f\n', ii, norm(abs(deltaf)), norm(u_bc));
         end
         status = false;
     end 
    
     % Update previous steps' meanings
     if (ii > 1) && (ii < imax)
        previousForce = deltaf;
        previousDisp = norm(u_bc);
     end
    
     

     