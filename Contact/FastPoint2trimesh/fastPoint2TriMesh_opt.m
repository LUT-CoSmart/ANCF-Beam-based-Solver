function [distances,project_pts,outside, projected_faces] = fastPoint2TriMesh_opt(inputs, pts)


    % --- Inputs / precompute ---
    faces = inputs.faces;
    nodes = inputs.nodes;

    need_norm = ~(isfield(inputs,'face_mean_nodes') && isfield(inputs,'face_normals'));
    if need_norm
        [face_mean_nodes, face_normals] = getFaceCenterAndNormals(faces, nodes);
    else
        face_mean_nodes = inputs.face_mean_nodes;
        face_normals     = inputs.face_normals;
    end

    if isfield(inputs,'tree_model')
        tree_model = inputs.tree_model;
    else
        % KDTreeSearcher is fine; createns is a bit faster on recent MATLABs.
        try
            tree_model = createns(face_mean_nodes,'NSMethod','kdtree'); % R2016b+
        catch
            tree_model = KDTreeSearcher(face_mean_nodes);               % fallback
        end
    end

    % --- Nearest face by in/centroid KD-tree ---
    near_id = knnsearch(tree_model, pts);
    projected_faces = near_id;

    % --- Gather triangle data for each point (Qx3) ---
    A = nodes(faces(near_id,1),:);
    B = nodes(faces(near_id,2),:);
    C = nodes(faces(near_id,3),:);
    n = face_normals(near_id,:);                 % per-point face normal (unit)

    % --- Project each point to the supporting plane of its nearest face ---
    % Use vertex A to define plane: proj = p - dot(p - A, n)*n
    Ap = pts - A;
    d_plane = sum(Ap .* n, 2);                  % signed distance to plane
    proj = pts - d_plane .* n;                  % initial projection points (Qx3)

    % --- Barycentric test (vectorized) for "inside triangle" ---
    % Solve proj = A + u*e1 + v*e2, with e1 = B-A, e2 = C-A
    e1 = B - A; e2 = C - A;
    v2 = proj - A;

    % Build 2x2 system in dot space: [dot11 dot12; dot12 dot22][u;v] = [dot1p; dot2p]
    dot11 = sum(e1.*e1, 2);
    dot22 = sum(e2.*e2, 2);
    dot12 = sum(e1.*e2, 2);
    dot1p = sum(e1.*v2, 2);
    dot2p = sum(e2.*v2, 2);

    denom = (dot11 .* dot22 - dot12 .* dot12);
    % Handle degenerate triangles (zero area): mark with denom==0
    degenerate = denom <= eps(max(1, denom));
    invDen = zeros(size(denom), 'like', denom);
    invDen(~degenerate) = 1 ./ denom(~degenerate);

    u = ( dot22 .* dot1p - dot12 .* dot2p) .* invDen;
    v = (-dot12 .* dot1p + dot11 .* dot2p) .* invDen;

    inside = ~degenerate & (u >= 0) & (v >= 0) & (u + v <= 1);

    % --- Edge projections for outside points (vectorized) ---
    % Helper to project P to segment X->Y: Q = X + clamp(t,0,1)*(Y-X), t = dot(P-X, Y-X)/||Y-X||^2
    function [Q, d2] = proj_to_seg(X, Y, P)
        V  = Y - X;
        VV = sum(V.*V, 2);
        % Avoid division by zero for zero-length edges:
        t  = zeros(size(VV), 'like', VV);
        nz = VV > 0;
        t(nz) = sum((P(nz,:) - X(nz,:)) .* V(nz,:), 2) ./ VV(nz);
        t = max(0, min(1, t));
        Q  = X + t .* V;
        D  = P - Q;
        d2 = sum(D.*D, 2);
    end

    if any(~inside)
        idx = ~inside;

        % Candidates: edges AB, BC, CA
        [Qab, d2ab] = proj_to_seg(A(idx,:), B(idx,:), proj(idx,:));
        [Qbc, d2bc] = proj_to_seg(B(idx,:), C(idx,:), proj(idx,:));
        [Qca, d2ca] = proj_to_seg(C(idx,:), A(idx,:), proj(idx,:));

        % Choose minimum-distance candidate per point
        % Compare three distances without loops:
        d2 = [d2ab, d2bc, d2ca];
        [~, k] = min(d2, [], 2);

        Qmin = Qab;                       % default
        sel  = (k == 2); Qmin(sel,:) = Qbc(sel,:);
        sel  = (k == 3); Qmin(sel,:) = Qca(sel,:);

        proj(idx,:) = Qmin;               % overwrite outside points with edge projection
    end

    % --- Signed distances: use face normal and original sign (outside/inside) ---
    % Original sign from plane distance (like your 'signs'):
    sgn = sign(d_plane);
    outside = sgn >= 0;

    % Geometric distances from pts to final projected points, keep sign
    distances = vecnorm(pts - proj, 2, 2) .* sgn;

    % --- Outputs ---
    project_pts = proj;
end
