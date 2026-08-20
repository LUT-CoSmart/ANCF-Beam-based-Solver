function visualizationSurfaceNormalStress(Body,Show)
    
    
    if Show 
        figure();
        axis equal
        visualization(Body,Body.q,'none',Show); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        xloc = Body.xloc;
        faces = Body.BodyFaces;
        IsoData = Body.IsoData;
        SurfacePoints = feval(Body.SurfacefunctionName, Body, Body.q); 
        
        getFaceCenterAndNormals = Body.getFaceCenterAndNormals;
        [~,face_normals]=getFaceCenterAndNormals(faces,SurfacePoints);

        points = [];

        for k = 1:size(faces,1)  
            
            normal = face_normals(k,:)';
            Xi_points = IsoData(faces(k,:),:);
            pos_points = SurfacePoints(faces(k,:),:);
            n = size(Xi_points, 1); % Determine the number of integration points (3 or 4)

            for i = 1:n
                r = pos_points(i,:)';
                xi = Xi_points(i,1);
                eta = Xi_points(i,2);
                zeta = Xi_points(i,3);
                Element = Xi_points(i,4);

                uk=Body.u(xloc(Element,:));  
                qk=Body.q(xloc(Element,:)); 
                qk0=Body.q0(xloc(Element,:)); 
                Phik=Body.Phim(Element,:)';

                F = Body.F(qk,qk0,uk,xi,eta,zeta);
                S = Body.S(F,qk0,Phik,xi,eta,zeta);

                Sigma = Body.Sigma(F,S); % Cauchy stress tensor 
                Sigma_surf = normal' * Sigma * normal;
                points = [points; r' Sigma_surf];
            end
        end

        h = scatter3(points(:,1), points(:,2), points(:,3),25,points(:,4),'filled','MarkerEdgeColor', 'k');

        h.DataTipTemplate.DataTipRows(4).Label = '\sigma_{nn}';
        h.DataTipTemplate.DataTipRows(4).Format = '%.6e';
        colormap(turbo);
        cb = colorbar;
        cb.Label.String = 'Surface stress';
        title('Surface stress');

    else
        disp('No visualization')
    end   