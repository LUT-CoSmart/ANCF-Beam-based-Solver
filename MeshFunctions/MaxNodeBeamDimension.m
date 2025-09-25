function R = MaxNodeBeamDimension(Body)

     dx = Body.Length.Ln/(2*(Body.ElementNodes-1)); %(  *--)(--*--)(--*  )
                                                    % 3-nodded element with closest distance for each node
                                                    %  half of the distance between nodes
     dyz = sqrt(Body.Length.Y^2 + Body.Length.Z^2)/2;
     R = max([dx, dyz],[],'all');
     % R = max([dx, Body.Length.Y, Body.Length.Z]);   