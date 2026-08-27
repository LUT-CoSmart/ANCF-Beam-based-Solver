function  [distances,project_pts,projected_faces,projected_normals,projected_elements,inside]=ProjectPoints(Body,Points)
           

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
           [distances,project_pts,outside,projected_faces]=Point2Mesh(inputs,Points); 
           projected_normals = -face_normals(projected_faces, :); % find the normals to the surfaces (chaning from inwards) 
           projected_elements = faceElem(projected_faces);  % find the elements for respected surfaces
           inside = ~outside;
