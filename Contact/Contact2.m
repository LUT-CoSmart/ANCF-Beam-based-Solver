function [Kc,Fc,Gap] = Contact2(Body1,Body2,ContactType,ContactVariable)

    if ContactType == "None"
       Fc = zeros(Body1.TotalDofs + Body2.TotalDofs,1);
       Kc = zeros(length(Fc));
       Gap = NaN;
    else        
       %% TODO: add boxing to identify the necessity of the contact, for now we always consider its existence
       h = 10^(-8);
       
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

               Body1.q = q1_backup - 2*h*I_vec(1:TotalDofs1); % -2h
               Fm2h = ContactForce(Body1, Body2, ContactVariable, ContactType);

               Body1.q = q1_backup -   h*I_vec(1:TotalDofs1);  % -h
               Fm1h = ContactForce(Body1, Body2, ContactVariable, ContactType);

               Body1.q = q1_backup +   h*I_vec(1:TotalDofs1);  % +h
               Fp1h = ContactForce(Body1, Body2, ContactVariable, ContactType);

               Body1.q = q1_backup + 2*h*I_vec(1:TotalDofs1);  % +2h
               Fp2h = ContactForce(Body1, Body2, ContactVariable, ContactType);                
           else                    
               Body2.q = q2_backup - 2*h*I_vec(1+TotalDofs1:TotalDofs);  
               Fm2h = ContactForce(Body1, Body2, ContactVariable, ContactType);

               Body2.q = q2_backup -   h*I_vec(1+TotalDofs1:TotalDofs);  
               Fm1h = ContactForce(Body1, Body2, ContactVariable, ContactType);
 
               Body2.q = q2_backup +    h*I_vec(1+TotalDofs1:TotalDofs);  
               Fp1h = ContactForce(Body1, Body2, ContactVariable, ContactType); 
                
               Body2.q = q2_backup + 2*h*I_vec(1+TotalDofs1:TotalDofs);  
               Fp2h = ContactForce(Body1, Body2, ContactVariable, ContactType);                
           end

           Kc(:,ii) = (-Fp2h + 8*Fp1h - 8*Fm1h + Fm2h) / (12*h);           
           % Kc(:,ii) = (Fm2h  - 4*Fm1h + 3*Fc) / (2*h);
           % Kc(:,ii) = (Fp1h - Fm1h) / (2*h);
           
           
           I_vec(ii)=0;    % returning all back to normal  
       end
       Body1.q = q1_backup; % restore 
       Body2.q = q2_backup; % restore


       % simple
       lambda = eps / h * norm(Kc, 'fro');   % version 2
       Kc = Kc + lambda * eye(TotalDofs);
    end       
