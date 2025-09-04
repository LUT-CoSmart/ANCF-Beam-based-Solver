function [Kc,Fc,Gap,GapMax] = Contact(Body1,Body2,ContactType,ContactVariable,ContactRegType)

    if ContactType == "None"
       Fc = zeros(Body1.TotalDofs + Body2.TotalDofs,1);
       Kc = zeros(length(Fc));
       Gap = NaN;
       GapMax.gap = 0;
       GapMax.area = NaN;
    else 
        % Creating names for the bodies' surface functions

        Body1.SurfacePoints = feval(Body1.SurfacefunctionName, Body1, Body1.q);         
        Body2.SurfacePoints = feval(Body2.SurfacefunctionName, Body2, Body2.q);
        %% TODO: add boxing to identify the necessity of the contact, for now we always consider its existence
        h = 2*sqrt(eps);
       
        TotalDofs1 = Body1.TotalDofs;
        TotalDofs2 = Body2.TotalDofs;
        TotalDofs = TotalDofs1 + TotalDofs2;
        

        % Initialize the global contact forces
        Kc = zeros(TotalDofs,TotalDofs);
        [Fc,Gap,GapMax] = ContactForce(Body1,Body2,ContactVariable,ContactType);
            
        % variation of the variables
        I_vec=zeros(TotalDofs,1);

        % Backup original coordinates
        q1_backup = Body1.q;
        q2_backup = Body2.q;
        SurfacePoints1_backup = Body1.SurfacePoints;
              
        u1_backup = Body1.u;
        u2_backup = Body2.u;
        SurfacePoints2_backup = Body2.SurfacePoints;

        for ii = 1:TotalDofs
            
           I_vec(ii)=1;

           % this split is to distribute coord. between bodies 
           if ii <= TotalDofs1

               Body1.q = q1_backup - h*I_vec(1:TotalDofs1); 
               Body1.u = u1_backup - h*I_vec(1:TotalDofs1);
               Body1.SurfacePoints = feval(Body1.SurfacefunctionName + "Partly", Body1, Body1.q, ii); 
               [Fch,~,~] = ContactForce(Body1,Body2,ContactVariable, ContactType); % force due to variation            
               Body1.SurfacePoints = SurfacePoints1_backup;
           else   
               % h = max(sqrtEps * abs(q2_backup(ii-TotalDofs1)) , h1); 
               
               Body2.q = q2_backup - h*I_vec(1+TotalDofs1:TotalDofs);  
               Body2.u = u2_backup - h*I_vec(1+TotalDofs1:TotalDofs);  
               Body2.SurfacePoints = feval(Body2.SurfacefunctionName + "Partly", Body2, Body2.q, ii - TotalDofs1); 
               [Fch, ~, ~] = ContactForce(Body1,Body2,ContactVariable, ContactType); % force due to variation            
               Body2.SurfacePoints = SurfacePoints2_backup; 
           end
          
           Kc(:,ii) = (Fc - Fch) / h;
           I_vec(ii)=0;   

        end

        Body1.q = q1_backup; % restore
        Body2.q = q2_backup; % restore

        Body1.u = u1_backup; % restore
        Body2.u = u2_backup; % restore
        
        [~,Kc] = Regularization(Kc,Fc,ContactRegType,false);
        
    end      
    
    