function [alphas,fs,gs] = Zoom(myFx,x0,d,alphal,alphah,fx0,gx0)
    % Algorithms 3.2 on page 59 in 
    % Numerical Optimization, by Nocedal and Wright
    % This function is called by strongwolfe
    c1 = 1e-4;
    c2 = 0.9;
    i =0;
    maxIter = 5;

    while true
       % bisection
       alphax = 0.5*(alphal+alphah);
       alphas = alphax;
       xx = x0 + alphax*d;
       [fxx,gxx] = feval(myFx,xx);
       fs = fxx;
       gs = gxx;
       gxx = gxx'*d;
       xl = x0 + alphal*d;
       fxl = feval(myFx,xl);
       if ((fxx > fx0 + c1*alphax*gx0) | (fxx >= fxl)),
          alphah = alphax;
       else
          if abs(gxx) <= -c2*gx0,
            alphas = alphax;
            return;
          end
          if gxx*(alphah-alphal) >= 0
            alphah = alphal;
          end
          alphal = alphax;
       end
         i = i+1;
       if i > maxIter
          alphas = alphax;
          return
       end
    end
end% end of Zoom