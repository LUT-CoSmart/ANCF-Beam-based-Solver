function Outcome = FindProjection(BodyToProject, Body)
        
        Outcome = [];
        fl = false; % initial contact possiblity, used later for function refactoring 
        tol = sqrt(eps);
        q = Body.q;
        Surfaces = Body.BodySurfaceType;
        Point2Mesh = Body.Point2Mesh;
        [PointsToProject,isoData,mask] = SelectingProjectedPoints_v2(BodyToProject, Body, tol);
        %  verification of direction (works for v1 without trisurf)        
        % [PointsToProject,isoData,mask,Nodes,facesToProject,centersSelected,normalsSelected,dirFace] = SelectingProjectedPoints_v2(BodyToProject, Body, tol);        
        % plot3(PointsToProject(:,1),PointsToProject(:,2),PointsToProject(:,3),'or' )
        % hold on
        % plot3(Body.SurfacePoints(:,1),Body.SurfacePoints(:,2),Body.SurfacePoints(:,3),'ok' )
        % quiver3(centersSelected(:,1), centersSelected(:,2), centersSelected(:,3),dirFace(:,1), dirFace(:,2), dirFace(:,3), 0.2, 'b', 'LineWidth', 1.5); 
        % plot3(Nodes(:,1),Nodes(:,2),Nodes(:,3),'-ok') 
        % plot3(BodyToProject.SurfacePoints(:,1),BodyToProject.SurfacePoints(:,2),BodyToProject.SurfacePoints(:,3),'og' )
        % trisurf(facesToProject, PointsToProject(:,1), PointsToProject(:,2), PointsToProject(:,3), 'FaceAlpha', 0.35,'EdgeColor', 'k');
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        

        if any(mask) % contact distantly possible
           fl = true; 
           getFaceCenterAndNormals = Body.getFaceCenterAndNormals;
           % Getting element parameters
           faces = Body.BodyFaces;
           faceElem = Body.BodyFacesElements;
           SurfacePoints = Body.SurfacePoints; 

           [face_mean_nodes,face_normals]=getFaceCenterAndNormals(faces,SurfacePoints);

           inputs.faces=faces;
           inputs.nodes=SurfacePoints;
           inputs.face_mean_nodes=face_mean_nodes;
           inputs.face_normals= -face_normals; % change normals for the distance calculation to outward,
                                               % because the calculated normals are inwards  

           [distances,~,outside,projected_faces]=Point2Mesh(inputs,PointsToProject); 
                      
           highlight_face = faces(projected_faces, :); % Get the vertex indices of the selected face  
           highlight_normals = face_normals(projected_faces, :); % find the normals to the surfaces
           idx = faceElem(projected_faces);  % find the elements for respected surfaces  

        end 
                  
        if fl
            Inside = ~outside;                               
            Inside(Inside) = abs(distances(Inside)) > tol;    
        else
            
        end    
        if (fl) && any(Inside) % choosing points inside, due to the normal identification procedure distances>0 
            
            xloc = Body.xloc;           
            Shape = Body.Shape;
            nabla_r_xi = Body.nabla_r_xi;
            distancesInside = distances(Inside);
            idxInside = idx(Inside);

            n = length(distancesInside);
            
            PointInside = PointsToProject(Inside,:);
            isoData = isoData(Inside,:);

            FaceNormal= highlight_normals(Inside,:);
            Face = highlight_face(Inside,:);
            
            % point isocoord targ (3) & element number (1) & distance (1) &
            % normal (3) & point isocoord cond (4 + element no. ) & contact area
            % contact point in global CS (3) 
            Outcome = zeros(n,16); % prelocation

            for i = 1:n
        
                qk=q(xloc(idxInside(i),:)); % current element number
          
                xi_eta_zeta_result = FindIsoCoord(Shape,nabla_r_xi,qk,PointInside(i,:));

                % patch area where the point is projected to (to calculate contact stresses?)                
                A =  SurfacePoints(Face(i,1),:)';
                B =  SurfacePoints(Face(i,2),:)';
                C =  SurfacePoints(Face(i,3),:)';

                if Surfaces == "triangles"
                    Area = 0.5 * norm( cross(B - A, C - A) );
                else
                    D =  SurfacePoints(Face(i,4),:)';
                    Area = 0.5 * norm(cross(B - D, C - A));
                end

                Outcome(i,:) =  [xi_eta_zeta_result', idxInside(i), distancesInside(i), FaceNormal(i,:), isoData(i,:), Area, PointInside(i,:)]; 

            end 
        end