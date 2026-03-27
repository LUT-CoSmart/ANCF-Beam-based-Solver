function [u_bc,deltaf] = Newton_Krylov(ii,Body,Fext,Re,name)
    % Taken from (7):
    % J. Hales, S. Novascone, R. L. Williamson, D. Gaston, M. Tonks (2012). 
    % Solving Nonlinear Solid Mechanics Problems with the Jacobian-Free Newton Krylov Method.
    % Computer Modeling in Engineering and Sciences. 84. 123. 
    
    if nargin < 5
        name = "JF";
        warning('chosen JF algorithm');  
    end

    persistent v r s Fe J;
           
    q_backup = Body.q;
    u_backup = Body.u;
    bc = Body.bc; 
    
    if ii == 1
        [J,Fe] = InnerForce(Body);
        r =  -(Fe - Fext);       % assembley    
        s = zeros(size(r));
        v = r;
        u_bc =  zeros(size(r(bc)));

    else
        u_bc =  zeros(size(r(bc)));
        nx = norm(q_backup(bc));        
        nv = norm(v(bc));
        h = 1e-7 * (1 + nx) / max(nv, 1e-20);
        
        r_bc = r(bc);
        v(~bc) = 0;
        v_bc = v(bc);
        
        switch name 

            case "JF"
                g = J*v;
    
            case "CG"
                Body.q = q_backup + h*v;
                Body.u = u_backup + h*v;
    
                [~,Fev] = InnerForce(Body);
                g = (Fev - Fe) / h;
    
                Body.q = q_backup; % restore
                Body.u = u_backup; % restore

            otherwise
                warning(' not correct Newton-Krylov algorithm, switched to JF');    
        end

        g_bc= g(bc);
        
        alpha = r_bc' * r_bc/(v_bc'*g_bc);
        s = s + alpha * v;
        r = r - alpha * g; % update r

        beta = r(bc)' * r(bc) / (r_bc' *r_bc);
        v = r + beta * v;
    end  

    deltaf=r(bc)/norm(Fext(Body.bc));
        
    if all(abs(deltaf) < Re) % to make the exit from the algorith similar to others
       u_bc = s(bc);
    end    