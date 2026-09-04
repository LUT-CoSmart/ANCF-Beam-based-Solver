function SurfaceShape = BuildBeamSurfaceMatrixLocal(Body)

    Shape_ = Body.Shape;
    xloc = Body.xloc;
    IsoData = Body.IsoData;
    N_points = size(IsoData,1);
    TotalDofs = Body.TotalDofs;
    ElementDofs = Body.ElementDofs;

    N_entries = Body.DIM * N_points * ElementDofs;

    I = zeros(N_entries,1);
    J = zeros(N_entries,1);
    V = zeros(N_entries,1);

    counter = 0;

    for i = 1:N_points

        xi    = IsoData(i,1);
        eta   = IsoData(i,2);
        zeta  = IsoData(i,3);
        Element = IsoData(i,4);
        DOFs  = xloc(Element,:);
        Shape = Shape_(xi,eta,zeta);

        rows = 3*(i-1) + (1:3); % rows corresponding to x,y,z of surface point

        for j = 1:3

            ind = counter + (1:ElementDofs);

            I(ind) = rows(j);
            J(ind) = DOFs;
            V(ind) = Shape(j,:);

            counter = counter + ElementDofs;

        end
    end

    SurfaceShape = sparse(I,J,V,Body.DIM*N_points,TotalDofs);


    