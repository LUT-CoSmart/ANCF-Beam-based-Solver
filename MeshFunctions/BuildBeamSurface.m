function SurfacePoints = BuildBeamSurface(Body,q)

        
        xloc = Body.xloc;
        Shape = Body.Shape;
        L = Body.Length.Ln;
        H = Body.Length.Y;
        W = Body.Length.Z;
        IsoData = Body.IsoData;
        n = size(IsoData,1);
        SurfacePoints = zeros(n,3);

        for i = 1:n
            xi = IsoData(i,1);
            eta = IsoData(i,2);
            zeta = IsoData(i,3);
            Element = IsoData(i,4);
            qk = q(xloc(Element ,:));
            r = Shape(L,H,W,xi,eta,zeta)*qk;
            SurfacePoints(i,:) = r';
        end
        

        


        
 
   


