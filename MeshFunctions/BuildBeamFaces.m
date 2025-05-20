function Body = BuildBeamFaces(Body)
         

        faces = [];
        
        switch Body.IntegrationType

           case {"Standard", "Poigen"} % fixed CS   
                
                pointCS = Body.SurfaceXi.pointCS;  
                Nxi = Body.SurfaceXi.Nxi;
                NpointCS = size(pointCS,1);

                % Creating faces with a triangle mesh      
                CSNumber = Body.ElementNumber * (Nxi-1) + 1; % total number of "Corss sections"

                for k = 1:CSNumber - 1
                    offset1 = (k-1)*NpointCS;
                    offset2 = k*NpointCS;
                    for p = 1:NpointCS - 1                
                        faces = [faces;
                                offset1+p,   offset1+p+1, offset2+p+1;
                                offset1+p,   offset2+p+1, offset2+p];
                    end
                    
                    % closing
                    faces = [faces;
                             offset1+NpointCS, offset1+1, offset2+1;
                             offset1+NpointCS, offset2+1, offset2+NpointCS];
                end    


                % side surfaces
                % All of it due to complex CS (like Tendon)
                edges = [(1:NpointCS)' [2:NpointCS 1]'];  % closed loop of edges
                dt = delaunayTriangulation(pointCS, edges); % creating nice triangulation to avoid collinearity for ending
                pgon = polyshape(pointCS); % polyshape from boundary points

                tri = dt.ConnectivityList;
                centroids = (pointCS(tri(:,1), :) + pointCS(tri(:,2), :) + pointCS(tri(:,3), :)) / 3;               
                in = isinterior(pgon, centroids(:,1), centroids(:,2)); % only triangles inside
                tri = tri(in, :);  % filtered triangles, leaving only inside the polygon


                startLayer = 1:NpointCS;
                endLayer   = (CSNumber - 1) * NpointCS + (1:NpointCS);

                facesStart = startLayer(tri);
                facesEnd   = endLayer(tri);

                facesStart = facesStart(:, [1 3 2]);  % inward at start
                facesEnd   = facesEnd(:, [1 2 3]);    % inward at end (adjust if needed)

                faces = [faces; facesStart; facesEnd];


            otherwise                  
                error('****** Unkown Integration type for %s ******', Body.Name);
        end
        
        Body.BodyFaces = faces;
        
        % For contact: in this way of the points' organization, the normals of the trimesh is directed to the inside volume 
        % Option        
        % SurfacePoints = BuildBeamSurface(Body,Body.q0);
        % addpath(genpath("Contact"))
        % [mean_nodes,face_normals]=getFaceCenterAndNormals(faces,SurfacePoints);
        % quiver3(mean_nodes(:,1), mean_nodes(:,2), mean_nodes(:,3),face_normals(:,1),  face_normals(:,2),  face_normals(:,3), 0.5, 'r', 'LineWidth', 0.01); 
        % patch('Vertices',SurfacePoints,'Faces',faces,'FaceColor','cyan','EdgeColor','black');