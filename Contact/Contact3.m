function [Kc,Fc,Gap] = Contact3(Body1,Body2,ContactType,ContactVariable)

    if ContactType == "None"
       Fc = zeros(Body1.TotalDofs + Body2.TotalDofs,1);
       Kc = zeros(length(Fc));
       Gap = NaN;
    else        
       %% TODO: add boxing to identify the necessity of the contact, for now we always consider its existence
       h = 10^(-9);
       h2 = h / 2;
       TotalDofs1 = Body1.TotalDofs;
       TotalDofs2 = Body2.TotalDofs;
       TotalDofs = TotalDofs1 + TotalDofs2;
       
       % Initialize the global contact forces
       Kc = zeros(TotalDofs,TotalDofs);

       [Fc,Gap] = ContactForce(Body1,Body2,ContactVariable,ContactType);
        

       
       
       % variation of the variables
       I_vec=zeros(TotalDofs,1);
       % Backup original coordinates
       q1_backup = Body1.q;
       q2_backup = Body2.q;
       for ii = 1:TotalDofs
           I_vec(ii)=1;

           % this split is to distribute coord. between bodies 
           if ii <= TotalDofs1

               Body1.q = q1_backup - h2*I_vec(1:TotalDofs1); % -2h
               Fmh2 = ContactForce(Body1, Body2, ContactVariable, ContactType);

               Body1.q = q1_backup -   h*I_vec(1:TotalDofs1);  % -h
               Fmh = ContactForce(Body1, Body2, ContactVariable, ContactType);

               Body1.q = q1_backup +   h*I_vec(1:TotalDofs1);  % +h
               Fph = ContactForce(Body1, Body2, ContactVariable, ContactType);

               Body1.q = q1_backup + h2*I_vec(1:TotalDofs1);  % +2h
               Fph2 = ContactForce(Body1, Body2, ContactVariable, ContactType); 

               Body1.q = q1_backup; % restore
           else                    
               Body2.q = q2_backup - h2*I_vec(1+TotalDofs1:TotalDofs);  
               Fmh2 = ContactForce(Body1, Body2, ContactVariable, ContactType);

               Body2.q = q2_backup -   h*I_vec(1+TotalDofs1:TotalDofs);  
               Fmh = ContactForce(Body1, Body2, ContactVariable, ContactType);
 
               Body2.q = q2_backup +   h*I_vec(1+TotalDofs1:TotalDofs);  
               Fph = ContactForce(Body1, Body2, ContactVariable, ContactType); 
                
               Body2.q = q2_backup + h2*I_vec(1+TotalDofs1:TotalDofs);  
               Fph2 = ContactForce(Body1, Body2, ContactVariable, ContactType); 

               Body2.q = q2_backup; % restore
           end

           D = (Fph - Fmh) / (2*h);
           D2 = (Fph2 - Fmh2) / (2*h2);

           Kc(:,ii) = (4 * D2- D) / 3;       
           % returning all back to normal
           I_vec(ii)=0;   
           
           
       end
        
    end       
