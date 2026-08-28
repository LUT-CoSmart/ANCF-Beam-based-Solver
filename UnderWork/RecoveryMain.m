function Recovery = RecoveryMain(Body)
        
        

        xloc = Body.xloc;
        TotalDofs = Body.TotalDofs;
        ElementDofs = Body.ElementDofs;
        ElementNumber = Body.ElementNumber;

        uk=Body.u(xloc(k,:));  
        qk = Body.q(xloc(k,:));
        qk0=Body.q0(xloc(k,:)); 
        Phik=Body.Phim(k,:)';

        Recovery=zeros(TotalDofs,1);
        for i = 1:ElementNumber
            Recovery_local = zeros(ElementDofs,1);
            Recovery_local =  Recovery_v2(Body,qk,q0k,uk,Phik);
            

            




            for jj = 1:ElementDofs
                ind01 = xloc(ii,jj); % global row
                Recovery(ind01) = Recovery(ind01) + Recovery_local(jj);
            end
        end    