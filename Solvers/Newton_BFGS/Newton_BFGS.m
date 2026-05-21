% function [u_bc,deltaf] = Newton_BFGS(step, Body, Fext)
% 
%     persistent Sm Ym 
% 
%     if step == 1
%        n = Body.nx; 
%        m = floor(n/2); 
%        Sm = zeros(n,m);
%        Ym = zeros(n,m);
% 
%        [K,Fe] = InnerForce(Body);
%        ff =  Fe - Fext;                 % assembley
%        ff_bc=ff(Body.bc);               % Eliminate linear constraints from force vector        
%        K_bc =  K(Body.bc,Body.bc);
%        [alpha,f1,g1] = strongwolfe(InnerForce,-K_bc,Body,ff_bc,K_bc);
%    end    
% 
% 
% end
% 
% %%%%%%%%%%%%%%%%%
% function [alphas,fs,gs] = strongwolfe(myFx,d,x0,fx0,gx0)
%     % Function strongwolfe performs Line search satisfying strong Wolfe conditions
%     % Input
%     %   myFx:   the optimized function handle
%     %   d:      the direction we want to search
%     %   x0:     vector of initial start
%     %   fx0:    the function value at x0
%     %   gx0:    the gradient value at x0  
%     % Output
%     %   alphas: step size
%     %   fs:     the function value at x0+alphas*d
%     %   gs:     the gradient value at x0+alphas*d
%     %
%     % Notice
%     %   I use f and g to save caculation. This funcion strongwolfe is called by LBFGS_opt.m.
%     % Ref
%     %   Numerical Optimization, by Nocedal and Wright
%     % Author:
%     %   Guipeng Li @THU
%     %   guipenglee@gmail.com
%     maxIter = 3;
%     alpham = 20;
%     alphap = 0;
%     c1 = 1e-4;
%     c2 = 0.9;
%     alphax = 1;
%     gx0 = gx0'*d;
%     fxp = fx0;
%     gxp = gx0;
%     i=1;
% 
%     % Line search algorithm satisfying strong Wolfe conditions
%     % Algorithms 3.2 on page 59 in Numerical Optimization, by Nocedal and Wright
%     % alphap is alpha_{i-1}
%     % alphax is alpha_i
%     % alphas is what we want.
%         while true
%           xx = x0 + alphax*d;
%           [fxx,gxx] = feval(myFx,xx);
%           fs = fxx;
%           gs = gxx;
%           gxx = gxx'*d;
%           if (fxx > fx0 + c1*alphax*gx0) | ((i > 1) & (fxx >= fxp))
%             [alphas,fs,gs] = Zoom(myFx,x0,d,alphap,alphax,fx0,gx0);
%             return;
%           end
%           if abs(gxx) <= -c2*gx0,
%             alphas = alphax;
%             return;
%           end
%           if gxx >= 0,
%             [alphas,fs,gs] = Zoom(myFx,x0,d,alphax,alphap,fx0,gx0);
%             return;
%           end
% 
%           alphap = alphax;
%           fxp = fxx;
%           gxp = gxx;
% 
%           if i > maxIter
%               alphas = alphax;
%               return
%           end  
%           % r = rand(1);%randomly choose alphax from interval (alphap,alpham)
%           r = 0.8;
%           alphax = alphax + (alpham-alphax)*r;
%           i = i+1;
%         end
% end% end of strongwolfe
% 
% %%%%%%%
% function [alphas,fs,gs] = Zoom(myFx,x0,d,alphal,alphah,fx0,gx0)
% % Algorithms 3.2 on page 59 in 
% % Numerical Optimization, by Nocedal and Wright
% % This function is called by strongwolfe
% c1 = 1e-4;
% c2 = 0.9;
% i =0;
% maxIter = 5;
% 
%     while true
%        % bisection
%        alphax = 0.5*(alphal+alphah);
%        alphas = alphax;
%        xx = x0 + alphax*d;
%        [fxx,gxx] = feval(myFx,xx);
%        fs = fxx;
%        gs = gxx;
%        gxx = gxx'*d;
%        xl = x0 + alphal*d;
%        fxl = feval(myFx,xl);
%        if ((fxx > fx0 + c1*alphax*gx0) | (fxx >= fxl)),
%           alphah = alphax;
%        else
%           if abs(gxx) <= -c2*gx0,
%             alphas = alphax;
%             return;
%           end
%           if gxx*(alphah-alphal) >= 0
%             alphah = alphal;
%           end
%           alphal = alphax;
%        end
%          i = i+1;
%        if i > maxIter
%           alphas = alphax;
%           return
%        end
%     end
% end% end of Zoom
% 
% %%%%%%%%%%%%%%%%
% function Hg = getHg_lbfgs(g,S,Y,hdiag)
% % This function returns the approximate inverse Hessian multiplied by the gradient, H*g
% % Input
% %   S:    Memory matrix (n by k) , s{i}=x{i+1}-x{i}
% %   Y:    Memory matrix (n by k) , df{i}=df{i+1}-df{i}
% %   g:    gradient (n by 1)
% %   hdiag value of initial Hessian diagonal elements (scalar)
% % Output
% %   Hg    the the approximate inverse Hessian multiplied by the gradient g
% % Notice
% % This funcion getHg_lbfgs is called by LBFGS_opt.m.
% % Ref
% %   Nocedal, J. (1980). "Updating Quasi-Newton Matrices with Limited Storage".
% %   Wiki http://en.wikipedia.org/wiki/Limited-memory_BFGS
% %   two loop recursion
% 
%     [n,k] = size(S);
%     for i = 1:k
%         ro(i,1) = 1/(Y(:,i)'*S(:,i));
%     end
% 
%     q = zeros(n,k+1);
%     r = zeros(n,1);
%     alpha =zeros(k,1);
%     beta =zeros(k,1);
% 
%     % step 1
%     q(:,k+1) = g;
% 
%     % first loop
%     for i = k:-1:1
%         alpha(i) = ro(i)*S(:,i)'*q(:,i+1);
%         q(:,i) = q(:,i+1)-alpha(i)*Y(:,i);
%     end
% 
%     % Multiply by Initial Hessian
%     r = hdiag*q(:,1);
% 
%     % second loop
%     for i = 1:k
%         beta(i) = ro(i)*Y(:,i)'*r;
%         r = r + S(:,i)*(alpha(i)-beta(i));
%     end
%     % 
%     Hg=r;
% end % end of getHg_lbfgs