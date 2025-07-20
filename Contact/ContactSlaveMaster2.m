function [Fcont, Ftarg, Gap] = ContactSlaveMaster2(ContactBody,TargetBody,ContactVariable, ContactType)
                                                      %(BodySlave,BodyMaster,ContactVariable, ContactType)                  

          functionName = "Build" + ContactBody.ElementType + "Surface"; 
          SurfacePointsSlave = feval(functionName, ContactBody, ContactBody.q);

          OutcomeOnBody = FindProjection(SurfacePointsSlave, ContactBody.IsoData, TargetBody);

          % 0. contact force & gap initialization  
          Gap = 0;
          Fcont = zeros(ContactBody.TotalDofs,1);
          Ftarg = zeros(TargetBody.TotalDofs,1);

          % exctrating information (1 & 2 are just to keep order, it doesn't necessary correlate with possible bodies' names)
          Shape_cont = ContactBody.Shape;  
          L_cont = ContactBody.Length.Ln;
          H_cont = ContactBody.Length.Y;
          W_cont = ContactBody.Length.Z;

          Shape_targ = TargetBody.Shape;  
          L_targ = TargetBody.Length.Ln;
          H_targ = TargetBody.Length.Y;
          W_targ = TargetBody.Length.Z;

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

                 gap = abs(Xi(i,5)); 
                 normal = Xi(i,6:8); % NB: it is not original normal, outwards respected to the body   

                 Normal_targ = -normal';
                 Normal_cont =  normal';
                
                 penalty = ContactVariable;                    
                 if ContactType == "Penalty"
                     
                    Fcont_loc =  penalty * gap * Normal_cont;                                                                              
                    Ftarg_loc =  penalty * gap * Normal_targ;

                 elseif ContactType == "NitscheLin"
                    
                    q_targ = TargetBody.q(DOFs_targ);
                    u_targ = TargetBody.u(DOFs_targ);
                    q0_targ = TargetBody.q0(DOFs_targ);
                    phi_targ=TargetBody.phim(Element_targ,:)';
                    q0PosDofs_targ = q0_targ(TargetBody.PosDofs);
                    F_targ = TargetBody.F(q_targ,u_targ,q0PosDofs_targ,phi_targ,L_targ,H_targ,W_targ,xi_targ,eta_targ,zeta_targ);        % Deformation gradient
                    Sigma_targ = TargetBody.Sigma(F_targ); 

                    q_cont = ContactBody.q(DOFs_cont);
                    u_cont = ContactBody.u(DOFs_cont);
                    q0_cont = ContactBody.q0(DOFs_cont);
                    phi_cont=ContactBody.phim(Element_cont,:)';
                    q0PosDofs_cont = q0_cont(ContactBody.PosDofs);

                    F_cont = ContactBody.F(q_cont,u_cont,q0PosDofs_cont,phi_cont,L_cont,H_cont,W_cont,xi_cont,eta_cont,zeta_cont);        % Deformation gradient
                    Sigma_cont = ContactBody.Sigma(F_cont); 
                    
                    Sigma_n = Normal_cont' * Sigma_cont * Normal_cont - Normal_targ' * Sigma_targ * Normal_targ;
                    lambda = gap^3 * norm(Sigma_n);

                    d_lambda_targ = norm(Sigma_n)*Normal_targ;
                    d_lambda_cont = norm(Sigma_n)*Normal_cont; 
                     
                    Ftarg_loc =  penalty * gap * Normal_targ + lambda * Normal_targ + 0 * gap * d_lambda_targ;
                    Fcont_loc =  penalty * gap * Normal_cont + lambda * Normal_cont + 0 * gap * d_lambda_cont;
                 else
                     error('****** Contact type is not implemneted ******')
                 end     
                    
                 Fcont(DOFs_cont) = Fcont(DOFs_cont) + Shape_cont(L_cont,H_cont,W_cont,xi_cont,eta_cont,zeta_cont)'*Fcont_loc;
                 Ftarg(DOFs_targ) = Ftarg(DOFs_targ) + Shape_targ(L_targ,H_targ,W_targ,xi_targ,eta_targ,zeta_targ)'*Ftarg_loc; 
             end                      
          end