function Body = CreateBC(Body, Force, Boundary) % Creates boundary conditions
    
    TotalDofs = Body.TotalDofs;
    DofsAtNode = Body.DofsAtNode;
    xlocName = 'xloc' + Body.ElementType;   
    DIM = Body.DIM;
    
    % Boundary Condition
    bc = true(1, TotalDofs); % Define vector of linear constraints  
    
    % Boundary position check
    if ~isfield(Boundary.Position, 'X')
       Boundary.Position.X = 0;
    end
    if ~isfield(Boundary.Position, 'Y')
       Boundary.Position.Y = 0;
    end
    if ~isfield(Boundary.Position, 'Z')
       Boundary.Position.Z = 0;
    end 

    NodalBC = FindNodalIDLocally(Body,Boundary.Position); 
    Dofs = 1:DofsAtNode; % taking all DOFs (in the beginning)

    switch Boundary.Type            
           case "full" 
                % all dofs are already considered

           case "reduced" % allowing contraction along the slopes 
                if Body.SubType == "ANCF" % only ANCF has slopes vectors to fix
                   if Body.Slope_x  
                      DofsOff = [DIM+1, 2*DIM + 2, 3*DIM + 3];
                   else
                     DofsOff = [DIM + 2, 2*DIM + 3];
                   end    
                   
                   Dofs = setdiff(Dofs, DofsOff); 
                else
                   error('****** This boundary condition is not recognized for %s ******', Body.SubType); 
                end   

           case "positions"
                Dofs = 1:DIM;

           case "none"
                Dofs = [];

           otherwise   
                error('****** The boundary condition type is not recognized ******');      
    end    

    bcInd = feval(xlocName,DofsAtNode,NodalBC,Dofs);
    if bcInd~=0
        bc(bcInd)=0;% number of degrees of freedom of system after linear constraints  
    end
    Body.bc = bc;
    Body.ndof = sum(bc);      % Number of unconstrained DOFs  
    
    % Force application    
    Fext = zeros(TotalDofs,1); % Initialize vector of ext forces
    NodalForce = FindNodalIDLocally(Body,Force.Position);
    fextInd = feval(xlocName,DofsAtNode,NodalForce,1:3);
      
    ForceVector = [Force.Maginutude.X; Force.Maginutude.Y; Force.Maginutude.Z];
    Fext(fextInd) = ForceVector;

    Body.Fext = Fext;
    Body.fextInd = NodalForce;
     