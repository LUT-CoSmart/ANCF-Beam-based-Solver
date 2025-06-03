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
        
        TotalDofs = Body.TotalDofs;
        xloc = Body.xloc;
        K=zeros(TotalDofs);
        Fint=zeros(TotalDofs,1); 
        ElementDofs = Body.ElementDofs;
        functionName = Body.ElementType + Body.SubType;  
        % for ii = 1:Body.ElementNumber               
        %     [K_loc,Fe] = feval(functionName,Body,ii);
        % 
        %     for jj = 1:ElementDofs
    	%         ind01 = xloc(ii,jj); %Index 01
        %         for kk = 1:ElementDofs
        %             ind02 = xloc(ii,kk); % Index 02
        %             K(ind01,ind02) = K(ind01,ind02)+K_loc(jj,kk);
        %         end
        %         Fint(ind01) = Fint(ind01)+Fe(jj); %%% sign infron of F
        %      end %End of Assembly pp
        % end
        
        ElementNumber = Body.ElementNumber;
        K_local_cell = cell(ElementNumber,1);
        F_local_cell = cell(ElementNumber,1);
        parfor ii = 1:ElementNumber      
            [K_loc, Fe] = feval(functionName, Body, ii);
            K_local_cell{ii} = K_loc;
            F_local_cell{ii} = Fe;
        end

        for ii = 1:ElementNumber
            K_loc = K_local_cell{ii};
            Fe = F_local_cell{ii};

            for jj = 1:ElementDofs
                ind01 = xloc(ii,jj); % global row
                for kk = 1:ElementDofs
                    ind02 = xloc(ii,kk); % global col
                    K(ind01,ind02) = K(ind01,ind02) + K_loc(jj,kk);
                end
                Fint(ind01) = Fint(ind01) + Fe(jj);
            end
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
        