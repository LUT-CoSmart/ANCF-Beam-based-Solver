function Fc = ContactForceWrapper(~, y, Body1, Body2, ContactVariable, ContactType)
    
    nq1 = Body1.TotalDofs;
    
    Body1.q = y(1:nq1);
    Body2.q = y(nq1+1:end);
    
    Body1.SurfacePoints = Body1.SurfacePointsFunction(Body1.q);
    Body2.SurfacePoints = Body2.SurfacePointsFunction(Body2.q);
    
    [Fc,~] = ContactForce(Body1,Body2,ContactVariable,ContactType);

    Fc = Fc(:);
end