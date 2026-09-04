function [Outcome,Inside] = FindProjection(BodyToProject, Body)
        
        Outcome = [];
        tol = 2*sqrt(eps);                
        [PointsToProject,isoData,mask,pointNormals] = SelectingProjectedPoints(BodyToProject, Body, tol);
        contactPossible = any(mask);

        if contactPossible % contact distantly possible                  
           [distances,project_pts,~,projected_normals,projected_el,Inside]=ProjectPoints(Body,PointsToProject);
        end 
        
        if any(Inside) % choosing points inside, due to the normal identification procedure distances>0 
            q = Body.q;
            xloc = Body.xloc;           
            Shape = Body.Shape;
            nabla_by_xi = Body.nabla_by_xi;
            distancesInside = distances(Inside);
            idxInside = projected_el(Inside); % suitable elements DOFs the face on which s projected 

            n = length(distancesInside);
            
            PointInside = PointsToProject(Inside,:);
            PointInsideNormals = pointNormals(Inside,:);
            isoData = isoData(Inside,:); % from BodyToProject
            project_pts = project_pts(Inside,:);
            
            Outcome = repmat(struct('PointProjected',[], ...
                                    'IsoCoorProjected',[], ...
                                    'ElementProjected',[], ...
                                    'ContactPoint',[],...
                                    'IsoCoordContactPoint',[], ...
                                    'ElementContactPoint',[], ...
                                    'Gap',[], ...
                                    'NormalProjected',[],...
                                    'ContactPointNormals',[]), n, 1); % prelocation

            FaceNormal = projected_normals(Inside,:);   
            % Surfaces = Body.BodySurfaceType;
            % highlight_face = faces(projected_faces, :); % Get the vertex indices of the selected face
            % Face = highlight_face(Inside,:);

            

            for i = 1:n        
                qk=q(xloc(idxInside(i),:)); % current element number
          
                isoCoordinates = FindIsoCoord(Shape,nabla_by_xi,qk,project_pts(i,:));

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
                Outcome(i).NormalProjected     = FaceNormal(i,:); 
                Outcome(i).ContactPointNormals = PointInsideNormals(i,:); 
            
            end
        end