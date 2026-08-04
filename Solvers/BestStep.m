function Body = BestStep(bestStep,lambda,Body,deltaf,ii,imax)
       
    if bestStep   

       persistent deltaf_best Body_best iteration 
       
       delta =  norm(abs(deltaf)); 
       
       if lambda == 1 && ii == 1 
       % lambda here is accounted only if bestStep activated and we are already on another step of backtracking
       % thus searching and taking the best step not within current lambda but all; 
       
          iteration = ii;
          deltaf_best = delta;
          Body_best = Body;

       else 

          if delta < deltaf_best
             deltaf_best = delta;
             Body_best = Body;
             iteration = ii;
          end

       end
    
       if ii==imax                
          Body = Body_best; 
          warning('Chosen the following solution: Iteration: %d, Convergence: %10.4f: %10.7f\n', iteration, deltaf_best); 
       end 

    end  