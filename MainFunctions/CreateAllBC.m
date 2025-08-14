function Body = CreateAllBC(Body, Force, Displacement, Boundary) % Creates boundary conditions
    
    TotalDofs = Body.TotalDofs;
    DofsAtNode = Body.DofsAtNode;
    xlocName = 'xloc' + Body.ElementType;   
    DIM = Body.DIM;
    
    % Boundary Condition
    %bc = true(1, TotalDofs); % Define vector of linear constraints  
    Applied_bc = zeros(1, TotalDofs);
          
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

    NodalFixed = FindNodalIDLocally(Body,Boundary.Position); 
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

    % Applied displacements on the non-fixed end

    NodalDispl = FindNodalIDLocally(Body,Displacement.Position); 

    DofsDispl = [];

    if isfield(Displacement.Maginutude, 'X')
        DofsDispl = [DofsDispl 1];
    end
    if isfield(Displacement.Maginutude, 'Y')
        DofsDispl = [DofsDispl 2];
    end
    if isfield(Displacement.Maginutude, 'Z')
        DofsDispl = [DofsDispl 3];
    end

    FixedInd = feval(xlocName,DofsAtNode,NodalFixed,Dofs);
    DisplInd = feval(xlocName,DofsAtNode,NodalDispl,DofsDispl);

    bcInd = [FixedInd DisplInd];


    if bcInd~=0
        bc(bcInd)=0;% number of degrees of freedom of system after linear constraints  
    end
    Body.bc = bc;
    Body.ndof = sum(bc);      % Number of unconstrained DOFs  
    
    % Applying non-zero displacements  
    
    DisplacementVector = [Displacement.Maginutude.X; Displacement.Maginutude.Y; Displacement.Maginutude.Z];
    Applied_bc(DisplInd) = DisplacementVector;
    
    Body.applied_disp = Applied_bc;
    Body.bcInd = bcInd;
    
    % Force application    
    Fext = zeros(TotalDofs,1); % Initialize vector of ext forces
    NodalForce = FindNodalIDLocally(Body,Force.Position);
    fextInd = feval(xlocName,DofsAtNode,NodalForce,1:3);
     
     % Force check
    if ~isfield(Force.Maginutude, 'X')
       Force.Maginutude.X = 0;
    end
    if ~isfield(Force.Maginutude, 'Y')
       Force.Maginutude.Y = 0;
    end
    if ~isfield(Force.Maginutude, 'Z')
       Force.Maginutude.Z = 0;
    end 
    
    ForceVector = [Force.Maginutude.X; Force.Maginutude.Y; Force.Maginutude.Z];
    
    % 
    Fext(fextInd) = ForceVector;
    Body.Fext = Fext;
    Body.fextInd = fextInd;
    
    Body.ForceVectorInit = ForceVector;    