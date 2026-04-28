clc,clear,close all;
format long

addpath(genpath(pwd));

Body.Name = "Body";

% ########### Problem data ################################################
Body = DefineElement(Body,"Beam","ANCF",3333,"None"); % 1 - BodyName, 2 - type (beam, plate, etc.), 3 - element name, 4 - modification name (None, EDG, etc.)  
                                                      % ANCF Beam: 3243, 3333, 3343, 3353, 3363, 34X3 (34103) 
Body = Geometry(Body,"ten_Sol_3","Standard", "Gauss");  % Cross Sections: Rectangular, Oval, C, Tendon, etc.
                                                      % Integration points of generating line: Gauss, Lobatto  
Body = Materials(Body,"GOH", "optimized_SEE");        % Material models: GOH (GOH, Amir), Neo-Hookean (Neo), 2- and 5- constant Mooney-Rivlin (Mooney2, Mooney5),  Kirhhoff-Saint-Venant (KS).
                                                      % Integration Scheme: Poigen, Standard
                                                     
% ########### Complicate geometry ######################§##################
% Shift
Body.Shift.X = 0;
Body.Shift.Y = 0;
Body.Shift.Z = 0;

% Rotation (in degrees)
Body.Rotation.X = 0;
Body.Rotation.Y = 0;
Body.Rotation.Z = 0;

% Twist
Body.Twist.initial_rot = 0;
Body.Twist.angle = 0; % in degrees
Body.Twist.ro = 0;

% ########## Create FE Model ##############################################
ElementNumber = 1;
Body = CreateFEM(Body,ElementNumber);

% ########## Calculation adjustments ######################################
Body.FiniteDiference= "AceGen"; % Calculation of FD: Matlab, Matlab_automatic, AceGen
Body.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body.DeformationType = "Finite"; % Deformation type: Finite, Small
Body = AddTensors(Body);

%% TODO: rebuild CreateMex, it addresses the wrong folder
Body.mex = false;

% ########## Boundary Conditions ##########################################
% Force 
Force.Maginutude.X = 2e6;  % Elongation

% Positioning applied locally to the Undefomred configuration
% Shift and curvature are accounted automaticaly)
Force.Position.X = Body.Length.X;  
Force.Position.Y = 0;  
Force.Position.Z = 0; 

% Boundaries (applied locally, shift and curvature are accounted automaticaly)
Boundary.Position.X = 0;  
Boundary.Position.Y = 0;
Boundary.Position.Z = 0;
Boundary.Type = "reduced"; % there are several types: full, reduced, positions, none

Body = CreateBC(Body, Force, Boundary); % Application of Boundary conditions
 
% % %####################### Solving ######################################## 
steps = 300;  % sub-loading steps
titertot=0;  
Re=10^(-4);           % Stopping criterion for residual
imax=20;              % Maximum number of iterations for Newton's method 

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

    Body = SaveResults(Body,i,"all"); % options: "all", "last", each by (number) 
end
% POST PROCESSING ###############################################
visDeformed = true;
visInitial = true;
PostProcessing(Body,visDeformed,visInitial);
% volume check
faces=Body.BodyFaces;
vertices_before = feval(Body.SurfacefunctionName, Body, Body.q0);
vertices_after  = feval(Body.SurfacefunctionName, Body, Body.q);

V_after = VolumeViaFaces(vertices_after, faces);
V_before = VolumeViaFaces(vertices_before, faces);

fprintf('Volume before: %10.12f; Volume after: %10.12f; Relative change: %10.12f \n', V_before, V_after, (V_after-V_before)/V_before)
CleanTemp(Body, true)