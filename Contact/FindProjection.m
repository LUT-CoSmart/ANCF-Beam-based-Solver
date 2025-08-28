function Outcome = FindProjection(PointsToProject, isoData, Body)

        Outcome = [];

        faces = Body.BodyFaces;
        SurfacePoints = feval("Build" + Body.ElementType + "Surface", Body, Body.q);
        SurfacePointsIso = Body.IsoData;

        % Findinding the closest nodes to a point 
        DofID = xlocBeam(Body.DofsAtNode,1:Body.NodeNumber,1:3);
        Nodes = reshape(Body.q(DofID), 3, []).'; % nodes positions ,reorginized to NodeNumberx3 
        [~, Distance] = knnsearch(Nodes, PointsToProject, 'K', 1);

        % Checking the contact possibility
        dyz = sqrt(Body.Length.Y^2 + Body.Length.Z^2)/2;
        RofBodyNodes = max([Body.Length.Ln/Body.ElementNodes, dyz],[],'all');
        
        cond = Distance < RofBodyNodes;
        PointsToProject = PointsToProject(cond,:);           
        isoData = isoData(cond,:);

        fl = false; % flag, showing the contact possiblity
        
        if any(cond)
           fl = true; % contact distantly possible

           % Getting element parameters
           q = Body.q;
           xloc = Body.xloc;           
           Shape = Body.Shape;
           ShapeXi = Body.ShapeXi;
           ShapeEta = Body.ShapeEta;
           ShapeZeta = Body.ShapeZeta;
          
           %% TODO: decrease the set of points (PointsToProject) by choosing only those, 
           %% which are on the same side with Nodes:  dir = Point' - NodalPoint; if dot(dir,r-NodalPoint)>0.

           [face_mean_nodes,face_normals]=getFaceCenterAndNormals(faces,SurfacePoints);
           inputs.faces=faces;
           inputs.nodes=SurfacePoints;
           inputs.face_mean_nodes=face_mean_nodes;
           inputs.face_normals= - face_normals; % change normal to outward
           [distances,~,~,projected_faces]=fastPoint2TriMesh(inputs,PointsToProject,0);         

           %% TODO: This is true for beam elements only 
           
           highlight_face = faces(projected_faces, :); % Get the vertex indices of the selected face  
           highlight_normals = face_normals(projected_faces, :); % find the normals to the surfaces
           idx = SurfacePointsIso(highlight_face,4); % now we need to find the element via finding the closest nodes, which define the element
        end 
          
        tol = 1e-6;
        Inside = (distances < 0) & (abs(distances) > tol);

        if (fl) && any(Inside) % choosing points inside, due to the normal identification procedure distances>0  
    
            distancesInside = distances(Inside);
            idxInside = idx(Inside);
            PointInside = PointsToProject(Inside,:);
            isoData = isoData(Inside,:);
            FaceNormal= highlight_normals(Inside,:);
            Face = highlight_face(Inside,:);
            xi_eta_zeta_Array = []; % point isocoord & element number & distance, such array, just to save much info as I want, later to make it zero array  

            for i = 1:length(distancesInside)
        
                qk=q(xloc(idxInside(i),:));
          
                xi_eta_zeta_result = FindIsoCoord(Shape,ShapeXi,ShapeEta,ShapeZeta,qk, PointInside(i,:)); % it will be used for Nitsche
                
                % patch area where the point is projected to
                % in the case, there are plans to calculate contact stresses
                A =  SurfacePoints(Face(i,1),:)';
                B =  SurfacePoints(Face(i,2),:)';
                C =  SurfacePoints(Face(i,3),:)';
                
                Area = 1/2 * norm( cross(B - A, C - A) );

                xi_eta_zeta_Array(i,:) =  [xi_eta_zeta_result', idxInside(i), distancesInside(i), FaceNormal(i,:), isoData(i,:), Area]; 
           end  
           Outcome = xi_eta_zeta_Array;
        end