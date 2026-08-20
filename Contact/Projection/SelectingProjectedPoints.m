function [PointsToProject,isoData,mask] = SelectingProjectedPoints(BodyToProject, Body, tol)
    
    PointsToProject = BodyToProject.SurfacePoints;
    isoData = BodyToProject.IsoData;
    facesToProject = BodyToProject.BodyFaces; 
    getFaceCenterAndNormals = BodyToProject.getFaceCenterAndNormals;

    % Findinding the closest node to a point 
    q = Body.q;    
    DofID = feval(Body.xlocFunction,Body.DofsAtNode,1:Body.NodeNumber,1:3);
    Nodes = reshape(q(DofID), 3, []).'; % nodes positions, reorginized to NodeNumberx3        
        
    [face_mean_projected,face_normals_projected]=getFaceCenterAndNormals(facesToProject,PointsToProject);

    % closest node to each face center
    [NodeInx, Distance] = knnsearch(Nodes, face_mean_projected, 'K', 1);
    ClosestNodesFace  = Nodes(NodeInx,:);

    face_normals_projected = -face_normals_projected;
    dirFace = ClosestNodesFace - face_mean_projected;
    
    % contact distance condition for faces
    cond = Distance < Body.NodeSphere; % if the faces are in correct direction, but far 
    cond2 = dot(dirFace, face_normals_projected, 2) > tol;

    mask = cond & cond2;
    facesSelected = facesToProject(mask,:);
    
    % reduced point set
    PointIds = unique(facesSelected(:), 'stable');    
    PointsToProject = PointsToProject(PointIds,:);
    isoData = isoData(PointIds,:);