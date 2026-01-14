function Fc = ContactForceWrapper(~, y, Body1, Body2, ContactVariable, ContactType)
    
    nq1 = Body1.TotalDofs;

    
    Body1.q = y(1:nq1);
    Body2.q = y(nq1+1:end);

    
    Body1.SurfacePoints = feval(Body1.SurfacefunctionName, Body1, Body1.q); 
    Body2.SurfacePoints = feval(Body2.SurfacefunctionName, Body2, Body2.q);
    
    [Fc,~,~] = ContactForce(Body1,Body2,ContactVariable,ContactType);

    % numjac expects a vector
    Fc = Fc(:);
end