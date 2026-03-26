function [u_bc,deltaf,step] = Newton_Krylov_CG(Body,Fext,imax,Re)
    % Taken from (7):
    % J. Hales, S. Novascone, R. L. Williamson, D. Gaston, M. Tonks (2012). 
    % Solving Nonlinear Solid Mechanics Problems with the Jacobian-Free Newton Krylov Method.
    % Computer Modeling in Engineering and Sciences. 84. 123. 
    
    [~,Fe] = InnerForce(Body);       
    q_backup = Body.q;
    u_backup = Body.u;

    bc = Body.bc;
    
    r =  -(Fe - Fext);       % assembley    
    s = zeros(size(r));
    v = r;

    nx = norm(q_backup(bc));
    nv = norm(v(bc));
    
    for step = 1:imax
        h = 1e-7 * (1 + nx) / max(nv, 1e-20);
        r_bc = r(bc);
        v(~bc) = 0;
        v_bc = v(bc);

        Body.q = q_backup + h*v;
        Body.u = u_backup + h*v;

        [~,Fev] = InnerForce(Body); 

        Body.q = q_backup; % restore
        Body.u = u_backup; % restore
        
        g = (Fev - Fe) / h;
        g_bc= g(bc);
        alpha = r_bc' * r_bc/(v_bc'*g_bc);
        s = s + alpha * v;
        r = r - alpha * g; % update r

        if norm(r(bc)) < Re * norm(Fext(bc))
           break 
        end

        beta = r(bc)' * r(bc) / (r_bc' *r_bc);
        v = r + beta * v;
        nv = norm(v(bc));
    end

    deltaf=r(bc)/norm(Fext(Body.bc));% Compute residua
    u_bc = s(bc);
    