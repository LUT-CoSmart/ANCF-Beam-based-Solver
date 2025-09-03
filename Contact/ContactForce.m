function [Fc, Gap, GapMax] = ContactForce(Body1,Body2,ContactVariable,ContactType)
        
        if (Body1.ContactRole == "slave") && (Body2.ContactRole == "master")
            
            [Fc1, Fc2, Gap, GapMax] = ContactSlaveMaster(Body1, Body2, ContactVariable, ContactType);

        elseif (Body1.ContactRole == "master") && (Body2.ContactRole == "slave")
            
            [Fc2, Fc1, Gap, GapMax] = ContactSlaveMaster(Body2, Body1, ContactVariable, ContactType);
           
        elseif (Body1.ContactRole == "master") && (Body2.ContactRole == "master")
            
            % Projection of Body2 on Body1
            [Fc2_1, Fc1_1, Gap1, Gap1Max] = ContactSlaveMaster(Body2, Body1, ContactVariable, ContactType);

            % Projection of Body1 on Body2
            [Fc1_2, Fc2_2, Gap2, Gap2Max] = ContactSlaveMaster(Body1, Body2, ContactVariable, ContactType);
                        
            Fc1 = Fc1_1 + Fc1_2;
            Fc2 = Fc2_1 + Fc2_2; 
            Gap = max(Gap1,Gap2);

            if Gap1Max.gap > Gap2Max.gap     
                GapMax = Gap1Max;
            else
                GapMax = Gap2Max;
            end    

        else

             error('****** Unkown case of contact roles ******')
        end     
      

        Fc = [Fc1; Fc2]; 
          