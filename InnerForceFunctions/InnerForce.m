function [K,Fint] = InnerForce(Body)
    
        % Accounting rotation 
        % NB: Element Inner Force Functions aren't derived for taking ridgid body rotation
        if isfield(Body, 'Rotation')              
            RotInv = Body.RotInv;
            Rot = Body.Rot;
            Body.q = RotInv * Body.q;
            Body.q0 = RotInv * Body.q0;
            Body.u = RotInv * Body.u;
            Body.q0f = RotInv * Body.q0f;
        end  
        
        xloc = Body.xloc;
        K=zeros(Body.TotalDofs,Body.TotalDofs);
        Fint=zeros(Body.TotalDofs,1); 
        ElementDofs = Body.ElementDofs;
                  
        for ii = 1:Body.ElementNumber      
            
            functionName = Body.ElementType + Body.SubType;  
            [K_loc,Fe] = feval(functionName,Body,ii);
    
            for jj = 1:ElementDofs
    	        ind01 = xloc(ii,jj); %Index 01
                for kk = 1:ElementDofs
                    ind02 = xloc(ii,kk); % Index 02
                    K(ind01,ind02) = K(ind01,ind02)+K_loc(jj,kk);
                end
                Fint(ind01) = Fint(ind01)+Fe(jj); %%% sign infron of F
             end %End of Assembly pp
        end

        % Returning back
        if isfield(Body, 'Rotation')
            Body.q = Rot* Body.q;
            Body.q0 = Rot* Body.q0;
            Body.u = Rot* Body.u;
            Body.q0f = Rot * Body.q0f;
            Fint = Rot * Fint;
            K = Rot * K * RotInv;
        end
        