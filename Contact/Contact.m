function [Kc,Fc,Gap] = Contact(Body1,Body2,ContactType,ContactVariable)

    if ContactType == "None"
       Fc = zeros(Body1.TotalDofs + Body2.TotalDofs,1);
       Kc = zeros(length(Fc));
       Gap = NaN;
    else        
       %% TODO: add boxing to identify the necessity of the contact, for now we always consider its existence
       h = 10^(-9);
       TotalDofs1 = Body1.TotalDofs;
       TotalDofs2 = Body2.TotalDofs;
       TotalDofs = TotalDofs1 + TotalDofs2;
        

       % Initialize the global contact forces
       Kc = zeros(TotalDofs,TotalDofs);
       [Fc,Gap] = ContactForce(Body1,Body2,ContactVariable,ContactType);
            
       % variation of the variables
       I_vec=zeros(TotalDofs,1);
       for ii = 1:TotalDofs
           I_vec(ii)=1;

           % this split is to distribute coord. between bodies 
           if ii <= TotalDofs1
               Body1.q = Body1.q - h*I_vec(1:TotalDofs1); 
           else                    
               Body2.q = Body2.q - h*I_vec(1+TotalDofs1:TotalDofs);        
           end

           [Fch, ~] = ContactForce(Body1,Body2,ContactVariable, ContactType); % force due to variation            
           Kc(:,ii) = (Fc - Fch) / h; 

           % returning all back to normal
           if ii <= TotalDofs1
              Body1.q = Body1.q + h*I_vec(1:TotalDofs1); 
           else                    
              Body2.q = Body2.q + h*I_vec(1+TotalDofs1:TotalDofs);
           end       

           I_vec(ii)=0;
       end
    end       
