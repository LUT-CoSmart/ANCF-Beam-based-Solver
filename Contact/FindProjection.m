function Outcome = FindProjection(PointsToProject, isoData, Body)

        Outcome = [];

        faces = Body.BodyFaces;
        faceElem = Body.BodyFacesElements;
        SurfacePoints = Body.SurfacePoints;
        SurfacePointsIso = Body.IsoData;

        % Findinding the closest nodes to a point 
        DofID = xlocBeam(Body.DofsAtNode,1:Body.NodeNumber,1:3);
        Nodes = reshape(Body.q(DofID), 3, []).'; % nodes positions ,reorginized to NodeNumberx3 
        [~, Distance] = knnsearch(Nodes, PointsToProject, 'K', 1);

        % Checking the contact possibility
        cond = Distance < Body.NodeSphere;
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
           % [distances,~,outside,projected_faces]=fastPoint2TriMesh(inputs,PointsToProject,0);         
           [distances,~,outside,projected_faces]=fastPoint2TriMesh_opt(inputs,PointsToProject); 
           
           %% TODO: This is true for beam elements only 
           % attributes of the target body (Body here)
           highlight_face = faces(projected_faces, :); % Get the vertex indices of the selected face  
           highlight_normals = face_normals(projected_faces, :); % find the normals to the surfaces

           idx = SurfacePointsIso(highlight_face,4); % now we need to find the element
           % idx =  faceElem(highlight_face);
        end 
          
        tol =sqrt(eps);
        Inside = ~outside;                               
        Inside(Inside) = abs(distances(Inside)) > tol;    
        
        if (fl) && any(Inside) % choosing points inside, due to the normal identification procedure distances>0  
    
            distancesInside = distances(Inside);
            idxInside = idx(Inside);
            PointInside = PointsToProject(Inside,:);
            isoData = isoData(Inside,:);
            FaceNormal= highlight_normals(Inside,:);
            Face = highlight_face(Inside,:);
            xi_eta_zeta_Array = []; % point isocoord & element number & distance, such array, just to save much info as I want, later to make it zero array  

            for i = 1:length(distancesInside)
        
                qk=q(xloc(idxInside(i),:)); % current element number
          
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