function visualization(Body,q,FaceColor,Show)
    
    
    if Show
        functionName = "Build" + Body.ElementType + "Surface"; 
        vertices = feval(functionName, Body, q);
        DofID = xlocBeam(Body.DofsAtNode,1:Body.NodeNumber,1:3);
        NodePositions = reshape(q(DofID), 3, []).'; % nodes positions
            
        hold on
        patch('Vertices', vertices, 'Faces', Body.BodyFaces, 'FaceColor', FaceColor, 'FaceAlpha', 0.3);
        plot3(vertices(:,1), vertices(:,2), vertices(:,3), '.k')    
        plot3(NodePositions(:,1), NodePositions(:,2), NodePositions(:,3), '-*k')
        view(3);  
    else
        disp('No visualization')
    end   


    