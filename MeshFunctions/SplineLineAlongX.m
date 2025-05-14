function [Ln,drdx] = SplineLineAlongX(xk,yk,zk,n)
    % Here we are to obtain the elements' length and slope_X vectors 
    pp_yk = csaps(xk, yk);     % Spline for y-coordinates
    pp_zk = csaps(xk, zk);     % Spline for z-coordinates                
    
    dpp_yk = fnder(pp_yk);     % Derivative of spline for y-coordinates
    dpp_zk = fnder(pp_zk);     % Derivative of spline for z-coordinates

    
    dy = fnval(dpp_yk, xk);
    dz = fnval(dpp_zk, xk);
    dx = ones(length(xk),1);
    
    drdx = [dx, dy(:), dz(:)];          % Tangents for the beam line
    drdx = drdx ./ vecnorm(drdx, 2, 2); % Normalize each tangent vector
                                 
    d = [xk' yk' zk'];                                            % Collection of the points
    % d = [xk fnval(pp_yk) fnval(pp_zk)]                          % Another way, with approxiamed values (they are different for a rough approx.)
    CS = cat(1,0,cumsum(sqrt(sum(diff(d,[],1).^2,2))));           % Cumulative sum of the distances between points
    Ln = CS(end) / n;                                             % Length of the elements
    
    
    

