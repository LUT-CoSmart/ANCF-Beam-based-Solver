function [Fcont, Ftarg, Gap] = ContactSlaveMaster(ContactBody,TargetBody,ContactVariable, ContactType)
                                                                     
          Gap.Total = 0;
          Gap.Area = NaN; % doesn't matter, when GapMax.gap == 0
          Gap.Maximum = 0;
          Gap.Points=[];  

          % exctrating information
          Shape_cont = ContactBody.Shape;  
          Shape_targ = TargetBody.Shape;  
         
          % Contact force & gap initialization  
          Fcont = zeros(ContactBody.TotalDofs,1);
          Ftarg = zeros(TargetBody.TotalDofs,1);
        
          %% TODO: adding the bounding boxing to check the contact at the first place  

          % ProjectionNitsche

          Outcome = FindProjection(ContactBody, TargetBody);
          

          % Checking the contact presence
          if ~isempty(Outcome)

             if isequal(ContactType, @Nitsche)                 
               Outcome  = GapGradient(ContactBody, TargetBody, Outcome);  
             end

             for i = 1:numel(Outcome)  % loop over all points
                 Data  = Outcome(i);   
                 CurrentGap = abs(Data.Gap);

                 [Fcont_loc, Ftarg_loc] = ContactType(ContactVariable, ContactBody, TargetBody, Data);
                 
                 Gap.Total = Gap.Total + CurrentGap;
                 if CurrentGap > Gap.Maximum
                      Gap.Maximum = CurrentGap;                      
                 end   
                 
                 Xi_cont = Data.IsoCoordContactPoint;
                 Xi_targ = Data.IsoCoorProjected;

                 DOFs_cont =  ContactBody.xloc(Data.ElementContactPoint,:); % associated DOFs   
                 DOFs_targ =  TargetBody.xloc(Data.ElementProjected,:);    
                                  
                 Fcont(DOFs_cont) = Fcont(DOFs_cont) + Shape_cont(Xi_cont(1),Xi_cont(2),Xi_cont(3))'*Fcont_loc;
                 Ftarg(DOFs_targ) = Ftarg(DOFs_targ) + Shape_targ(Xi_targ(1),Xi_targ(2),Xi_targ(3))'*Ftarg_loc; 
                    
                 Gap.Points = [Gap.Points; Data.PointProjected Data.ElementProjected];   
             end                 
          end