function Outcome = GapGradient(BodyToProject, Body, Outcome)

n = numel(Outcome); % number of points of interest
h = 1e-7; % it is applied in local isoparametric coordinates 
          % and variates alongside with element size (no need being large)    

% 
Gaps = vertcat(Outcome.Gap);

% points to variate (from Body to Project)
q = BodyToProject.q;
xloc = BodyToProject.xloc;
Shape_ = BodyToProject.Shape;

IsoCoordContactPoints = vertcat(Outcome.IsoCoordContactPoint);
Element = vertcat(Outcome.ElementContactPoint);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Point2Mesh = Body.Point2Mesh;
getFaceCenterAndNormals = Body.getFaceCenterAndNormals;
faces = Body.BodyFaces;
SurfacePoints = Body.SurfacePoints; 
[face_mean_nodes,face_normals]=getFaceCenterAndNormals(faces,SurfacePoints);
inputs.faces=faces;
inputs.nodes=SurfacePoints;
inputs.face_mean_nodes=face_mean_nodes;
inputs.face_normals= -face_normals; 

ContactPoints =    vertcat(Outcome.ContactPoint)
[distances,outside,projected_faces]=Point2Mesh(inputs,ContactPoints);

Outcome2 = FindProjection(BodyToProject, Body);

vertcat(Outcome2.ContactPoint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

ContactPoints_variation = zeros(n,3);
for i = 1:3
    h_vec = -h*ones(n,1);
    IsoCoordContactPoints_variation = IsoCoordContactPoints; % Reset all coordinates 
    IsoCoordContactPoints_variation(:,i) = IsoCoordContactPoints(:,i) + h_vec;

    for j = 1:n % making sure the variation of xi is within the same element
        if IsoCoordContactPoints_variation(j,i) <-1
           h_vec(j) = -h_vec(j); 
           IsoCoordContactPoints_variation(j,i) = IsoCoordContactPoints(j,i) + h_vec(j); 
        end  
        
        xi = IsoCoordContactPoints_variation(j,1);
        eta = IsoCoordContactPoints_variation(j,2);
        zeta = IsoCoordContactPoints_variation(j,3);
        qk = q(xloc(Element(j) ,:));
        
        ContactPoints_variation(j,:) = qk' * Shape_(xi,eta,zeta)';
    end    
    

   


    [distances,outside,projected_faces]=Point2Mesh(inputs,ContactPoints)
    Normals_variation = -face_normals(projected_faces(~outside), :);
    
end    

    