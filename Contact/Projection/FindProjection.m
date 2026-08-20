function Outcome = FindProjection(BodyToProject, Body)
        
        Outcome = [];
        tol = 2*sqrt(eps);
        q = Body.q;
        % Surfaces = Body.BodySurfaceType;
        [PointsToProject,isoData,mask] = SelectingProjectedPoints(BodyToProject, Body, tol);

        contactPossible = any(mask);

        if contactPossible % contact distantly possible
           Point2Mesh = Body.Point2Mesh;
           getFaceCenterAndNormals = Body.getFaceCenterAndNormals;           
           faces = Body.BodyFaces; % Getting element parameters
           faceElem = Body.BodyFacesElements;
           SurfacePoints = Body.SurfacePoints; 

           [face_mean_nodes,face_normals]=getFaceCenterAndNormals(faces,SurfacePoints);

           inputs.faces=faces;
           inputs.nodes=SurfacePoints;
           inputs.face_mean_nodes=face_mean_nodes;
           inputs.face_normals= -face_normals; % change normals for the distance calculation to outward,
                                               % because the calculated normals are inwards  

           [distances,project_pts,outside,projected_faces]=Point2Mesh(inputs,PointsToProject); 
                      
           %highlight_face = faces(projected_faces, :); % Get the vertex indices of the selected face  
           highlight_normals = face_normals(projected_faces, :); % find the normals to the surfaces (inwards)
           idx = faceElem(projected_faces);  % find the elements for respected surfaces

           Inside = ~outside;
        end 
                  
        if any(Inside) % choosing points inside, due to the normal identification procedure distances>0 
            
            xloc = Body.xloc;           
            Shape = Body.Shape;
            nabla_r_xi = Body.nabla_r_xi;
            distancesInside = distances(Inside);
            idxInside = idx(Inside); % suitable elements DOFs the face on which s projected 

            n = length(distancesInside);
            
            PointInside = PointsToProject(Inside,:);
            isoData = isoData(Inside,:); % from BodyToProject

            FaceNormal= -highlight_normals(Inside,:); % chaning the sign of normal 
            % Face = highlight_face(Inside,:);
            project_pts = project_pts(Inside,:);

            Outcome = repmat(struct('PointProjected',[], ...
                                    'IsoCoorProjected',[], ...
                                    'ElementProjected',[], ...
                                    'ContactPoint',[],...
                                    'IsoCoordContactPoint',[], ...
                                    'ElementContactPoint',[], ...
                                    'Gap',[], ...
                                    'Normal',[]), n, 1); % prelocation

            for i = 1:n        
                qk=q(xloc(idxInside(i),:)); % current element number
          
                isoCoordinates = FindIsoCoord(Shape,nabla_r_xi,qk,project_pts(i,:));

                % patch area where the point is projected to (to calculate contact stresses?)                
                % A =  SurfacePoints(Face(i,1),:)';
                % B =  SurfacePoints(Face(i,2),:)';
                % C =  SurfacePoints(Face(i,3),:)';
                % if Surfaces == "triangles"
                %     Area = 0.5 * norm( cross(B - A, C - A) );
                % else
                %     D =  SurfacePoints(Face(i,4),:)';
                %     Area = 0.5 * norm(cross(B - D, C - A));
                % end

                % Store the results
                Outcome(i).PointProjected      = project_pts(i,:); 
                Outcome(i).IsoCoorProjected    = isoCoordinates;
                Outcome(i).ElementProjected    = idxInside(i); 
                Outcome(i).ContactPoint        = PointInside(i,:);
                Outcome(i).IsoCoordContactPoint= isoData(i,1:3);
                Outcome(i).ElementContactPoint = isoData(i,4);                               
                Outcome(i).Gap                 = distancesInside(i);
                Outcome(i).Normal              = FaceNormal(i,:)'; 
            end 
        end