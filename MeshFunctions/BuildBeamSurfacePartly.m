function SurfacePoints = BuildBeamSurfacePartly(Body,q,j)
        
        Shape_ = Body.Shape;
        xloc = Body.xloc;
        IsoData = Body.IsoData;
        SurfacePoints = Body.SurfacePoints; % before the change    
        
        % j - a changed position; 

        AffectedElements = find(any(xloc == j, 2)); % find DOFs with number j, it can be several elements 
        idx_points = find(ismember(IsoData(:,4), AffectedElements)); % all points belonging to these elements
        for i = 1:length(idx_points)            
            xi   = IsoData(idx_points(i),1);
            eta  = IsoData(idx_points(i),2);
            zeta = IsoData(idx_points(i),3);
            elems = IsoData(idx_points(i),4);
            q_affected   = q(xloc(elems,:));  % DOFs of the element 
            r = Shape_(xi,eta,zeta)*q_affected;
            SurfacePoints(idx_points(i),:) = r'; % update position
        end    
    
        


        
 
   


