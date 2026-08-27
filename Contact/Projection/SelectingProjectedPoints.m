function [PointsToProject,isoData,mask,pointNormals] = SelectingProjectedPoints(BodyToProject, Body, tol)

    PointsToProject = BodyToProject.SurfacePoints;
    isoData = BodyToProject.IsoData;
    facesToProject = BodyToProject.BodyFaces;
    getFaceCenterAndNormals = BodyToProject.getFaceCenterAndNormals;

    % Find the closest node to a point
    q = Body.q;

    DofID = feval(Body.xlocFunction, Body.DofsAtNode, 1:Body.NodeNumber,1:3);

    Nodes = reshape(q(DofID),3,[]).';

    [face_mean_projected,face_normals_projected] = getFaceCenterAndNormals(facesToProject,PointsToProject);

    % Closest node to each face center
    [NodeInx,Distance] = knnsearch( ...
        Nodes,face_mean_projected,'K',1);

    ClosestNodesFace = Nodes(NodeInx,:);

    face_normals_projected = -face_normals_projected;
    dirFace = ClosestNodesFace - face_mean_projected;

    % Contact distance condition for faces
    cond = Distance < Body.NodeSphere;
    cond2 = dot(dirFace,face_normals_projected,2) > tol;

    mask = cond & cond2;

    facesSelected = facesToProject(mask,:);
    selectedFaceNormals = face_normals_projected(mask,:);

    % Create this BEFORE reducing PointsToProject
    pointNormalsAll = zeros(size(PointsToProject,1),3);

    % Add each face normal to the points of that face
    for i = 1:size(facesSelected,1)
        ids = facesSelected(i,:);
        pointNormalsAll(ids,:) = pointNormalsAll(ids,:) + selectedFaceNormals(i,:);
    end

    % Reduced point set
    PointIds = unique(facesSelected(:),'stable');
    PointsToProject = PointsToProject(PointIds,:);
    isoData = isoData(PointIds,:);
    pointNormals = pointNormalsAll(PointIds,:);

    % Normalize point normals
    normalLength = vecnorm(pointNormals,2,2);
    pointNormals = pointNormals ./ max(normalLength,eps);