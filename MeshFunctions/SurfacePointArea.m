function SurfacePointArea = SurfacePointArea(Body, q)

    Faces    = Body.BodyFaces;
    N_points = size(Body.IsoData,1);
    Surface = Body.SurfacePointsFunction(q);   
    SurfacePointArea = zeros(N_points,1);

    NfaceNodes = size(Faces,2);

    for ff = 1:size(Faces,1)

        ids = Faces(ff,:);

        if NfaceNodes == 3
            
                x1 = Surface(ids(1),:);
                x2 = Surface(ids(2),:);
                x3 = Surface(ids(3),:);

                Aface = 0.5 * norm(cross(x2-x1, x3-x1));

                SurfacePointArea(ids) = SurfacePointArea(ids) + Aface/3;

        elseif NfaceNodes == 4

                x1 = Surface(ids(1),:);
                x2 = Surface(ids(2),:);
                x3 = Surface(ids(3),:);
                x4 = Surface(ids(4),:);

                A1 = 0.5 * norm(cross(x2-x1, x3-x1));
                A2 = 0.5 * norm(cross(x3-x1, x4-x1));

                Aface = A1 + A2;

                SurfacePointArea(ids) = SurfacePointArea(ids) + Aface/4;
        else
                error('Unsupported surface face');
        end
    end
