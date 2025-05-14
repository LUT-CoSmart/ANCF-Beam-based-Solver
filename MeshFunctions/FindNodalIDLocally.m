function NodalID = FindNodalIDLocally(Body,Position)
    
    
    tol=2*sqrt(eps); % calculation error

    P = Body.P00; % straight and undeformed configuration

    % check locations 
    if isfield(Position, 'X'), x = Position.X; else, x = 'all'; end
    if isfield(Position, 'Y'), y = Position.Y; else, y = 'all'; end
    if isfield(Position, 'Z'), z = Position.Z; else, z = 'all'; end


    mask = true(size(P,1),1); % autamotically accounted for "all"

    if isnumeric(x)
        mask = mask & abs(P(:,1) - x) < tol;
    end
    if isnumeric(y)
        mask = mask & abs(P(:,2) - y) < tol;
    end
    if isnumeric(z)
        mask = mask & abs(P(:,3) - z) < tol;
    end

    % Now mask contains only valid points
    NodalID = find(mask);
    
