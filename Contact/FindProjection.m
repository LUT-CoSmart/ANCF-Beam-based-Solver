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

           L = Body.Length.Ln;
           H = Body.Length.Y;
           W = Body.Length.Z;  

         
           %% TODO: decrease the set of points (PointsToProject) choosing only those, 
           %% which are on the same side with Nodes:  dir = Point' - NodalPoint; if dot(dir,r-NodalPoint)>0.

           [face_mean_nodes,face_normals]=getFaceCenterAndNormals(faces,SurfacePoints);
           inputs.faces=faces;
           inputs.nodes=SurfacePoints;
           inputs.face_mean_nodes=face_mean_nodes;
           inputs.face_normals= - face_normals; % change normal to outward
           [distances,~,~,projected_faces]=fastPoint2TriMesh(inputs,PointsToProject,1);         

           %% TODO: This is true for beam elements only 
           
           highlight_face = faces(projected_faces, :); % Get the vertex indices of the selected face  
           highlight_normals = face_normals(projected_faces, :); % find the normals to the surfaces
           idx = SurfacePointsIso(highlight_face(:, 1),4); % now we need to find the element via finding the closest nodes, which define the element
        end 
                
        if (fl) && any(distances<0) % choosing points inside, due to the normal identification procedure distances>0  
    
            Inside = distances<0;
            
            distancesInside = distances(Inside);
            idxInside = idx(Inside);
            PointInside = PointsToProject(Inside,:);
            isoData = isoData(Inside,:);
            FaceNormal= highlight_normals(Inside,:);

            xi_eta_zeta_Array = []; % point isocoord & element number & distance, such array, just to save much info as I want, later to make it zero array  

            for i = 1:length(distancesInside)
        
                qk=q(xloc(idxInside(i),:));

                % xi_eta_zeta0 = [0; 0; 0];                                        
                % options = optimoptions('lsqnonlin', 'Display','off','FunctionTolerance', 1e-5);
                % residual = @(xi_eta_zeta) Shape(L,H,W,xi_eta_zeta(1),xi_eta_zeta(2),xi_eta_zeta(3)) * qk - PointInside(i,:)'; % Defining the residual function
                % xi_eta_zeta_result = lsqnonlin(residual, xi_eta_zeta0, [-1; -1; -1], [1; 1; 1],options);
                % 
                xi_eta_zeta_result = FindIsoCoord(Shape,ShapeXi,ShapeEta,ShapeZeta,L,H,W,qk, PointInside(i,:));

                xi_eta_zeta_Array(i,:) =  [xi_eta_zeta_result', idxInside(i), distancesInside(i), FaceNormal(i,:), isoData(i,:)]; 
           end  
           Outcome = xi_eta_zeta_Array;
        end