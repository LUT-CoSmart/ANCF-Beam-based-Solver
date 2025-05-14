clc,clear,close all;
format long
addpath("MainFunctions");
addpath("Postprocessing");
addpath('MeshFunctions');
Body.Name = "Body";
% ########### Problem data ################################################
Body = DefineElement(Body,"Beam","ANCF",3333,"None");  % 1 - BodyName, 2 - type (beam, plate, etc.), 3 - element name, 4 - modification name (None, EDG, etc.)  
                                                       % ANCF Beam: 3243, 3333, 3343, 3353, 3363, 34X3 (34103)    
Body = Materials(Body,"GOH"); % Material models: Gas.-Ogd.-Hol. (GOH), Neo-Hookean (Neo), 2- and 5- constant Mooney-Rivlin (Mooney2, Mooney5),  Kirhhoff-Saint-Venant (KS).
% Itegration Scheme: Poigen, Standard
Body = Geometry(Body,"Rectangular","Standard");  % Cross Sections: Rectangular, Oval, C, Tendon
% ########### Complicate geometry #########################################
% Shift
Body.Shift.X = 0;
Body.Shift.Y = 0;
Body.Shift.Z = 0;

% Rotation (in degrees)
Body.Rotation.X = 0;
Body.Rotation.Y = 0;
Body.Rotation.Z = 0;

% Twist
Body.Twist.angle = 0; % in degrees
Body.Twist.ro = 0;
% ########## Create FE Model ##############################################
ElementNumber = 1;
Body = CreateFEM(Body,ElementNumber);
% ########## Calculation adjustments ######################################
Body.FiniteDiference= "Matlab"; % Calculation of FD: Matlab, AceGen
Body.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body.DeformationType = "Finite"; % Deformation type: Finite, Small
Body = AddTensors(Body);
% ########## Boundary Conditions ##########################################
% Force 
Force.Maginutude.X = 1e8;  % Elongation
Force.Maginutude.Y = 0;  
Force.Maginutude.Z = 0;  

% Positioning applied locally to the Undefomred configuration
% Shift and curvature are accounted automaticaly)
Force.Position.X = Body.Length.X;  
Force.Position.Y = 0;  
Force.Position.Z = 0; 

% Boundaries (applied locally, shift and curvature are accounted automaticaly)
Boundary.Position.X = 0;  
Boundary.Position.Y = 0;
Boundary.Position.Z = 0;
Boundary.Type = "full"; % there are several types: full, reduced, positions, none

% ########## Visualization of initial situation ###########################
Results = [];  
% % %####################### Solving ######################################## 
steps = 20;  % sub-loading steps
titertot=0;  
Re=10^(-4);                   % Stopping criterion for residual
imax=20;                      % Maximum number of iterations for Newton's method 

%START NEWTON'S METHOD   
for i=1:steps

    % Update forces, supported loading types: linear, exponential, quadratic, cubic;
    Subforce = SubLoading(Force, i, steps, "exponential"); 

    % Application of Boundary conditions
    Body = CreateBC(Body, Subforce, Boundary);
    Fext = Body.Fext;
      
              
    for ii=1:imax    
        tic; 
                                
        [K,Fe] = InnerForce(Body);
                       
        K_bc = K(Body.bc,Body.bc);            % Eliminate linear constraints from stiffness matrix
        ff =  Fe - Fext;

        ff_bc=ff(Body.bc);               % Eliminate linear constraints from force vector
        deltaf=ff_bc/norm(Fext(Body.bc));% Compute residual
        u_bc = -K_bc\ff_bc;             % Compute displacements
        Body.u(Body.bc) = Body.u(Body.bc)+u_bc;         % Add displacement to previous one
        Body.q(Body.bc) = Body.q(Body.bc)+u_bc;         % change the global positions
        titer=toc;
        titertot=titertot+titer;   
        if  all(abs(deltaf) < Re) || (norm(u_bc)<Re^2) 
            fprintf('Convergence: %10.4f, Displacements norm: %10.4f\n', norm(abs(deltaf)), norm(u_bc));
            fprintf('Solution is found on %d iteration, Total CPU-time: %f\n', ii, titertot);            
            break
        elseif ii==imax 
            fprintf('The solution is not found. The maximum number of iterations is reached. Total CPU-time: %d\n', ii);
        else     
            fprintf('Iteration: %d, Convergence: %10.4f, Displacements norm: %10.5f\n', ii, norm(abs(deltaf)), norm(u_bc));
        end              
    end           

      

    %Pick nodal displacements from result vector
    xlocName = 'xloc' + Body.ElementType;
    DofID = feval(xlocName,Body.DofsAtNode,Body.fextInd,1:3);
    uf = Body.u(DofID); 

    Results = [Results; Body.ElementNumber Body.TotalDofs uf'];
end
% POST PROCESSING ###############################################
visDeformed = true;
visInitial = true;
PostProcessing(Body,Results,visDeformed,visInitial) 