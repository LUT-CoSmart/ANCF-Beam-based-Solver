function Body = CreateFEM(Body,ElementNumber)
    
    
    Body.ElementNumber = ElementNumber;
    
    DofsAtNode = Body.DofsAtNode;
    DIM = Body.DIM;
    PosDofs = []; % Identification of positional DoFs within the chosen element 
    for i = 1:Body.ElementNodes
        if Body.Slope_x
            PosDofs = [PosDofs (i-1)*DofsAtNode+1:(i-1)*DofsAtNode+2*DIM];
        else   
            PosDofs = [PosDofs (i-1)*DofsAtNode+1:(i-1)*DofsAtNode+DIM];
        end    
    end   
    Body.PosDofs = PosDofs;

    MeshingFunctionName = Body.ElementType + Body.SubType + 'Mesh';
    Body = feval(MeshingFunctionName,Body);

    % gathering physical and geometrical element data in one vector
    Body.Dvec=[Body.const,Body.Length.Y,Body.Length.Z,Body.Length.Ln];
           
    % creating surface data for contact and visualization
    functionName = "Surface" + Body.ElementType + "Approx"; 
    Body = feval(functionName,Body);

    % creating faces and iso data of the whole body
    functionName = "Build" + Body.ElementType + "Faces";     
    Body = feval(functionName,Body); 
  
    
    
    
    