function [PointsToProject,isoData,mask,Nodes,facesToProject,centersSelected,normalsSelected,dirFace] = SelectingProjectedPoints_v2(BodyToProject, Body, tol)
    
    PointsToProject = BodyToProject.SurfacePoints;
    isoData = BodyToProject.IsoData;
    facesToProject = BodyToProject.BodyFaces; 

    % Findinding the closest node to a point 
    q = Body.q;    
    DofID = xlocBeam(Body.DofsAtNode,1:Body.NodeNumber,1:3); 
    Nodes = reshape(q(DofID), 3, []).'; % nodes positions, reorginized to NodeNumberx3        

    [face_mean_projected,face_normals_projected]=getFaceCenterAndNormals(facesToProject,PointsToProject);

    % closest target node to each face center
    [NodeInx, Distance] = knnsearch(Nodes, face_mean_projected, 'K', 1);
    ClosestNodesFace  = Nodes(NodeInx,:);

    % contact distance condition for faces
    cond = Distance < Body.NodeSphere;

    face_normals_projected = -face_normals_projected;
    dirFace = ClosestNodesFace - face_mean_projected;

    cond2 = dot(dirFace, face_normals_projected, 2) > tol;

    % mask = cond & cond2;
    mask = cond2;   
    
    facesSelected = facesToProject(mask,:);
    
    % reduced point set
    PointIds = unique(facesSelected(:), 'stable');    
    PointsToProject = PointsToProject(PointIds,:);
    isoData = isoData(PointIds,:);

    [~, facesToProject] = ismember(facesSelected, PointIds);    
    dirFace = dirFace(mask,:);
    centersSelected = face_mean_projected(mask,:);
    normalsSelected = face_normals_projected(mask,:);    