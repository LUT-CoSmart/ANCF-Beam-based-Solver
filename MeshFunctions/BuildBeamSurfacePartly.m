function SurfacePoints = BuildBeamSurfacePartly(Body,q,j)
        
        Shape_ = Body.Shape;
        xloc = Body.xloc;
        IsoData = Body.IsoData;
        SurfacePoints = Body.SurfacePoints; % before the change    
        
        % j - a changed position; 
        AffectedElements = find(any(xloc == j, 2));
        
        for i = 1:length(AffectedElements)
            
            AffectedElement = AffectedElements(i);
            AffectedDOFs =  xloc(AffectedElement,:); % affeted DOFs
            q_affected = q(AffectedDOFs);
            idx_affected_points = find(IsoData(:,4) == AffectedElement);

            for k = 1:length(idx_affected_points)

                point_iso = idx_affected_points(k);
                xi = IsoData(point_iso,1);
                eta = IsoData(point_iso,2);
                zeta = IsoData(point_iso,3);
                r = Shape_(xi,eta,zeta)*q_affected;
                SurfacePoints(point_iso,:) = r'; % update position

            end    
        end    
        


           
             
            


        


        
 
   


