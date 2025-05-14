function [FcSlave, FcMaster, Gap] = ContactSlaveMaster(BodySlave,BodyMaster,ContactVariable, ContactType)

          functionName = "Build" + BodySlave.ElementType + "Surface"; 
          SurfacePointsSlave = feval(functionName, BodySlave, BodySlave.q);

          OutcomeOnBody = FindProjection(SurfacePointsSlave, BodySlave.IsoData, BodyMaster);

          % 0. contact force & gap initialization  
          Gap = 0;
          FcSlave = zeros(BodySlave.TotalDofs,1);
          FcMaster = zeros(BodyMaster.TotalDofs,1);

          % exctrating information (1 & 2 are just to keep order, it doesn't necessary correlate with possible bodies' names)
          Shape1 = BodySlave.Shape;  
          L1 = BodySlave.Length.Ln;
          H1 = BodySlave.Length.Y;
          W1 = BodySlave.Length.Z;

          Shape2 = BodyMaster.Shape;  
          L2 = BodyMaster.Length.Ln;
          H2 = BodyMaster.Length.Y;
          W2 = BodyMaster.Length.Z;

          % Checking the contact presence
          if ~isempty(OutcomeOnBody)

             Xi = OutcomeOnBody;
             Gap = abs(sum(Xi(:,5)));

             for i = 1:size(Xi,1)
                 % forces to body1 from points under the SurfacePoints 
                 xi1 = Xi(i,9);
                 eta1 = Xi(i,10);
                 zeta1 = Xi(i,11);
                 Element1 = Xi(i,12);       % element of Body2 
                 DOFs1 =  BodySlave.xloc(Element1,:);     % associated DOFs
 
                 % forces to body2 from points of SurfacePoints inside
                 xi2 = Xi(i,1);
                 eta2 = Xi(i,2);
                 zeta2 = Xi(i,3);

                 Element2 = Xi(i,4);  % element  
                 DOFs2 =  BodyMaster.xloc(Element2,:);     % associated DOFs

                 gap = abs(Xi(i,5)); 
                 normal = Xi(i,6:8); % NB: it is not original normal, outwards respected to the body   

                 if ContactType == "Penalty"
                    penalty = ContactVariable;
                    FcSlave(DOFs1) = FcSlave(DOFs1) + penalty*gap*Shape1(L1,H1,W1,xi1,eta1,zeta1)'*normal';
                    FcMaster(DOFs2) = FcMaster(DOFs2) - penalty*gap*Shape2(L2,H2,W2,xi2,eta2,zeta2)'*normal'; 
                 else
                     error('****** Contact type is not implemneted ******')
                 end                     
             end                      
          end