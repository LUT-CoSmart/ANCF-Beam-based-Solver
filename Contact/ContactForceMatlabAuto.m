function [Kc,Fc,Gap,GapMax] = ContactForceMatlabAuto(Body1,Body2,ContactType,ContactVariable,CalculateStiffness)
    
    [Fc,Gap,GapMax] = ContactForce(Body1,Body2,ContactVariable,ContactType);
    
    Kc = zeros(length(Fc));

    if CalculateStiffness   
        q1 = Body1.q(:);
        q2 = Body2.q(:);
        y0 = [q1; q2];
        G = @(t,y) ContactForceWrapper(t,y,Body1,Body2,ContactVariable,ContactType);
        fac = 1e-4;
        Kc = numjac(G, 0, y0, Fc, fac, []);
    end