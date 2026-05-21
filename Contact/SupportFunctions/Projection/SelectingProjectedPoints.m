function [PointsToProject,isoData,PointNormals,mask,Nodes,dirPoint] = SelectingProjectedPoints(BodyToProject, Body, tol)
    
    PointsToProject = BodyToProject.SurfacePoints;

    isoData = BodyToProject.IsoData;
    facesToProject = BodyToProject.BodyFaces; 

    % Findinding the closest node to a point 
    q = Body.q;    
    DofID = xlocBeam(Body.DofsAtNode,1:Body.NodeNumber,1:3); 
    Nodes = reshape(q(DofID), 3, []).'; % nodes positions, reorginized to NodeNumberx3        
    [NodeInx, Distance] = knnsearch(Nodes, PointsToProject, 'K', 1); % taking the closest one 
    ClosestNodes = Nodes(NodeInx,:);
    
    cond = Distance < Body.NodeSphere;

    % Finding Normal for each point to project (see below why)
    [~,face_normals_projected]=getFaceCenterAndNormals(facesToProject,PointsToProject);
    [~, nVertFace] = size(facesToProject);
    ids = reshape(facesToProject.', [], 1);
    
    PointsToProject = PointsToProject(ids,:);
    isoData = isoData(ids,:);
    ClosestNodes = ClosestNodes(ids,:);

    PointNormals = -repelem(face_normals_projected, nVertFace, 1);

    % Decrease the set of points (PointsToProject) by choosing only those, 
    % which are on the same side with Nodes:  dir = Point' - NodalPoint; if dot(dir,r-NodalPoint)>0.       
    dirPoint = ClosestNodes - PointsToProject; % Direction from closest node to point
    cond2 = dot(dirPoint, PointNormals, 2) > tol;

    % mask = cond & cond2;
    mask = cond2;

    PointsToProject = PointsToProject(mask,:);
    isoData = isoData(mask,:);    
    PointNormals = PointNormals(mask,:);
    dirPoint = dirPoint(mask,:);
     
    [~, ia] = unique(PointsToProject, 'rows', 'stable');    % remove repeated points    
    PointsToProject = PointsToProject(ia,:);
    isoData = isoData(ia,:); 
    PointNormals = PointNormals(ia,:);
    dirPoint = dirPoint(ia,:);