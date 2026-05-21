function Hg = getHg_lbfgs(g,S,Y,hdiag)
    % This function returns the approximate inverse Hessian multiplied by the gradient, H*g
    % Input
    %   S:    Memory matrix (n by k) , s{i}=x{i+1}-x{i}
    %   Y:    Memory matrix (n by k) , df{i}=df{i+1}-df{i}
    %   g:    gradient (n by 1)
    %   hdiag value of initial Hessian diagonal elements (scalar)
    % Output
    %   Hg    the the approximate inverse Hessian multiplied by the gradient g
    % Notice
    % This funcion getHg_lbfgs is called by LBFGS_opt.m.
    % Ref
    %   Nocedal, J. (1980). "Updating Quasi-Newton Matrices with Limited Storage".
    %   Wiki http://en.wikipedia.org/wiki/Limited-memory_BFGS
    %   two loop recursion

    [n,k] = size(S);
    for i = 1:k
        ro(i,1) = 1/(Y(:,i)'*S(:,i));
    end

    q = zeros(n,k+1);
    r = zeros(n,1);
    alpha =zeros(k,1);
    beta =zeros(k,1);

    % step 1
    q(:,k+1) = g;

    % first loop
    for i = k:-1:1
        alpha(i) = ro(i)*S(:,i)'*q(:,i+1);
        q(:,i) = q(:,i+1)-alpha(i)*Y(:,i);
    end

    % Multiply by Initial Hessian
    r = hdiag*q(:,1);

    % second loop
    for i = 1:k
        beta(i) = ro(i)*Y(:,i)'*r;
        r = r + S(:,i)*(alpha(i)-beta(i));
    end
    % 
    Hg=r;
end % end of getHg_lbfgs