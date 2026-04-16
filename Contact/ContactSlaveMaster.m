function [Fcont, Ftarg, Gap, GapMax] = ContactSlaveMaster(ContactBody,TargetBody,ContactVariable, ContactType)
                                                                     
          GapMax.gap = 0;
          GapMax.area = NaN; % doesn't matter, when GapMax.gap == 0

          % exctrating information
          Shape_cont = ContactBody.Shape;  
          Shape_targ = TargetBody.Shape;  
         
          % Contact force & gap initialization  
          Gap = 0; % total gap
          Fcont = zeros(ContactBody.TotalDofs,1);
          Ftarg = zeros(TargetBody.TotalDofs,1);
        
          %% TODO: adding the bounding boxing to check the contact at the first place  

          % Projection
          Outcome = FindProjection(ContactBody.SurfacePoints, ContactBody.IsoData, TargetBody);
                  
          % Checking the contact presence
          if ~isempty(Outcome)

             for i = 1:size(Outcome,1)  % loop over all points
                                                     
                 [Fcont_loc, Ftarg_loc, DOFs_cont, DOFs_targ, Xi_cont, Xi_targ, gap] = ContactType(ContactVariable, ContactBody, TargetBody, Outcome(i,:));
                 
                 Gap = Gap + gap;
                 
                 if gap > GapMax.gap
                      GapMax.gap = gap;
                      GapMax.area = Outcome(i,13); 
                 end   
                 
                 % cross weighting                 
                 total_weight = norm(Fcont_loc) + norm(Ftarg_loc); 
                 weight_targ = norm(Fcont_loc) / total_weight;
                 weight_cont = norm(Ftarg_loc) / total_weight;
                 
                 Fcont(DOFs_cont) = Fcont(DOFs_cont) + Shape_cont(Xi_cont(1),Xi_cont(2),Xi_cont(3))'*Fcont_loc * weight_cont;
                 Ftarg(DOFs_targ) = Ftarg(DOFs_targ) + Shape_targ(Xi_targ(1),Xi_targ(2),Xi_targ(3))'*Ftarg_loc * weight_targ; 

             end                      
          end