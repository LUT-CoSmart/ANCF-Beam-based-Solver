clc,clear,close all;
format long
addpath("MainFunctions");
addpath("Postprocessing");
addpath('MeshFunctions');
addpath(genpath("Solvers"))
Body.Name = "Body";
% ########### Problem data ################################################
Body = DefineElement(Body,"Beam","ANCF",3333,"None");  % 1 - BodyName, 2 - type (beam, plate, etc.), 3 - element name, 4 - modification name (None, EDG, etc.)  
                                                       % ANCF Beam: 3243, 3333, 3343, 3353, 3363, 34X3 (34103)
% Itegration Scheme: Poigen, Standard
Body = Geometry(Body,'Rectangular',"Standard", "Gauss");  % Cross Sections: Rectangular, Oval, C, Tendon                                                       
                                                      % Integration points of generating line: Gauss, Lobatto      
Body = Materials(Body,'Neo'); % Material models: Gas.-Ogd.-Hol. (GOH), Neo-Hookean (Neo), 2- and 5- constant Mooney-Rivlin (Mooney2, Mooney5),  Kirhhoff-Saint-Venant (KS).

% ########### Complicate geometry #########################################
% Shift
Body.Shift.X = 0;
Body.Shift.Y = 0;
Body.Shift.Z = 0;

Body.Twist.angle = 0;
Body.Twist.initial_rot = 0;
Body1.Twist.ro = 0;
% Rotation (in degrees)
Body.Rotation.X = 0;
Body.Rotation.Y = 0;
Body.Rotation.Z = 0;

% Twist
Body.Twist.angle = 45; % in degrees
Body.Twist.ro = 0;

% ########## Create FE Model ##############################################
ElementNumber = 2;
Body = CreateFEM(Body,ElementNumber);

% ########## Calculation adjustments ######################################
Body.FiniteDiference= "AceGen"; % Calculation of FD: Matlab, AceGen
Body.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body.DeformationType = "Finite"; % Deformation type: Finite, Small
Body = AddTensors(Body);

%% TODO: rebuild CreateMex, it addresses the wrong folder
Body.mex = false;
% ########## Boundary Conditions ##########################################
% Force 
Force.Maginutude.X = 1e5;  % Elongation
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

Body = CreateBC(Body, Force, Boundary); % Application of Boundary conditions

% ########## Visualization of initial situation ########################### 
visualization(Body,Body.q0,'red',false);

% ####################### Solving ######################################## 
steps = 20;  % sub-loading steps
titertot=0;  
Re=10^(-5);                   % Stopping criterion for residual
imax=20;                      % Maximum number of iterations for Newton's method 

%START NEWTON'S METHOD   
for i=1:steps

    % Update forces, supported loading types: linear, exponential, quadratic, cubic;
    Body = SubLoading(Body, i, steps, "cubic"); 

    Fext = Body.Fext;    
    for ii=1:imax    
        tic; 
          
        [u_bc,deltaf] = Newton_full(Body,Fext);

        Body.u(Body.bc) = Body.u(Body.bc)+u_bc;         % Add displacement to previous one
        Body.q(Body.bc) = Body.q(Body.bc)+u_bc;         % change the global positions
        titer=toc;
        titertot=titertot+titer;   
        
        if printStatus(deltaf, u_bc, Re, i, ii, imax, steps, titertot)
            break;  
        end 

    end           

    Body = SaveResults(Body,i, "last"); % options: "all", "last", each by (number) 
end
% POST PROCESSING ###############################################
visDeformed = true;
visInitial = true;
PostProcessing(Body,visDeformed,visInitial) 
CleanTemp(Body, true)