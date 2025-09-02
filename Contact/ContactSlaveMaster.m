function [Fcont, Ftarg, Gap, GapMax] = ContactSlaveMaster(ContactBody,TargetBody,ContactVariable, ContactType)
                                                      %(BodySlave,BodyMaster,ContactVariable, ContactType)                  
          GapMax.gap = 0;
          GapMax.area = NaN; % doesn't matter, when GapMax.gap == 0

          OutcomeOnBody = FindProjection(ContactBody.SurfacePoints, ContactBody.IsoData, TargetBody);

          % 0. contact force & gap initialization  
          Gap = 0;
          Fcont = zeros(ContactBody.TotalDofs,1);
          Ftarg = zeros(TargetBody.TotalDofs,1);

          % exctrating information (1 & 2 are just to keep order, it doesn't necessary correlate with possible bodies' names)
          Shape_cont = ContactBody.Shape;  
          Shape_targ = TargetBody.Shape;  

          % Checking the contact presence
          if ~isempty(OutcomeOnBody)

             Xi = OutcomeOnBody;       % taking data 
             Gap = abs(sum(Xi(:,5)));  % total gap

             for i = 1:size(Xi,1)  % loop over all points

                 %  Data of slave (contact) body points under the surface of master body 
                 xi_cont = Xi(i,9);
                 eta_cont = Xi(i,10);
                 zeta_cont = Xi(i,11);
                 Element_cont = Xi(i,12);       % element of master body 
                 DOFs_cont =  ContactBody.xloc(Element_cont,:);     % associated DOFs
 
                 % Data of master (target) body points projected from slave ones
                 xi_targ = Xi(i,1);
                 eta_targ = Xi(i,2);
                 zeta_targ = Xi(i,3);
                 Element_targ = Xi(i,4);  % element  
                 DOFs_targ =  TargetBody.xloc(Element_targ,:);     % associated DOFs

                 Gap = abs(Xi(i,5)); 
                 
                 if Gap > GapMax.gap
                      GapMax.gap = Gap;
                      GapMax.area = Xi(i,13); 
                 end

                 Normal = Xi(i,6:8)'; % NB: it is not original normal, outwards respected to the body, so direction towards Contact   

                 Normal_targ = -Normal;
                 Normal_cont =  Normal;                

                 penalty = ContactVariable;                    
                 if ContactType == "Penalty"
                     
                    Fcont_loc =  penalty * Gap * Normal_cont;                                                                              
                    Ftarg_loc =  penalty * Gap * Normal_targ;

                 elseif ContactType == "NitscheLin"
                          
                    q_targ = TargetBody.q(DOFs_targ);
                    u_targ = TargetBody.u(DOFs_targ);
                    q0_targ = TargetBody.q0(DOFs_targ);
                    phi_targ=TargetBody.phim(Element_targ,:)';
                    q0PosDofs_targ = q0_targ(TargetBody.PosDofs);
                    F_targ = TargetBody.F(q_targ,u_targ,q0PosDofs_targ,phi_targ,xi_targ,eta_targ,zeta_targ);        % Deformation gradient
                    Sigma_targ = TargetBody.Sigma(F_targ); 

                    q_cont = ContactBody.q(DOFs_cont);
                    u_cont = ContactBody.u(DOFs_cont);
                    q0_cont = ContactBody.q0(DOFs_cont);
                    phi_cont=ContactBody.phim(Element_cont,:)';
                    q0PosDofs_cont = q0_cont(ContactBody.PosDofs);
                    F_cont = ContactBody.F(q_cont,u_cont,q0PosDofs_cont,phi_cont,xi_cont,eta_cont,zeta_cont);        % Deformation gradient
                    Sigma_cont = ContactBody.Sigma(F_cont); 
                                     
                    % Normal force difference 
                    Sigma_n = Normal_cont' * Sigma_cont * Normal_cont - Normal_targ' * Sigma_targ * Normal_targ;    
                    Gap_power = Gap;       
                    lambda = Gap_power * norm(Sigma_n);
                   
                    d_lambda_targ = norm(Sigma_n)*Normal_targ;
                    d_lambda_cont = norm(Sigma_n)*Normal_cont; 

                    Ftarg_loc = penalty * Gap * Normal_targ + lambda * Normal_targ + Gap_power * d_lambda_targ;
                    Fcont_loc = penalty * Gap * Normal_cont + lambda * Normal_cont + Gap_power * d_lambda_cont; 
                    
                 else
                     error('****** Contact type is not implemneted ******')
                 end     
                    
                 Fcont(DOFs_cont) = Fcont(DOFs_cont) + Shape_cont(xi_cont,eta_cont,zeta_cont)'*Fcont_loc;
                 Ftarg(DOFs_targ) = Ftarg(DOFs_targ) + Shape_targ(xi_targ,eta_targ,zeta_targ)'*Ftarg_loc; 

             end                      
          end