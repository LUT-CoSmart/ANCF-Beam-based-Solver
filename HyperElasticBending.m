clc,clear,close all;
format long
addpath(genpath(pwd));
Body.Name = "Body";
% ########### Problem data ################################################
Body = DefineElement(Body,"Beam","ANCF",3333,"None");  % 1 - BodyName, 2 - type (beam, plate, etc.), 3 - element name, 4 - modification name (None, EDG, etc.)  
                                                       % ANCF Beam: 3243, 3333, 3343, 3353, 3363, 34X3 (34103)    
% ANCF Beam: 3243, 3333, 3343, 3353, 3363, 34X3 (34103)
Body = DefineElement(Body,"Beam","ANCF",3363,"None");  
% Geometry
Body = Geometry(Body,"Rectangular","Standard","Gaus");  % Cross Sections: Rectangular, Oval, C, Tendon
                                                      % Integration points of generating line : Gauss, Lobatto
% Material models: GOH (GOH), Neo-Hookean (Neo), 2- and 5- constant Mooney-Rivlin (Mooney2, Mooney5),  Kirhhoff-Saint-Venant (KS).
Body = Materials(Body,'Neo',"Sol_old"); 
% ########## Create FE Model ##############################################
ElementNumber = 32;
Body = CreateFEM(Body,ElementNumber);

% % ########## Calculation adjustments ######################################
Body.FiniteDiference= "AceGen"; % Calculation of FD: Matlab, Matlab_automatic, AceGen
Body.SolutionBase = "Displacement"; % Solution-based calculation: Position, Displacement
Body.DeformationType = "Finite"; % Deformation type: Finite, Small

Body = AddTensors(Body);

% ########## Boundary Conditions ##########################################
% Force (applied locally, shift and curvature are accounted automaticaly)
Force.Maginutude.Y = -1e5;
Force.Position.X = Body.Length.X;  % Elongation

% Boundaries (applied locally, shift and curvature are accounted automaticaly)
Boundary.Position = [];
Boundary.Type = "full"; % there are several types: full, reduced, positions, none


% %####################### Solving ######################################## 
steps = 50;  % sub-loading steps, a lot for non full Newton-based algorithms
titertot=0; 

Body = CreateBC(Body, Force, Boundary); % Application of Boundary conditions

%START NEWTON'S METHOD
for i=1:steps

    % Update forces
    Body = SubLoading(Body, i, steps, "linear"); 

    Re=10^(-5);                   % Stopping criterion for residual
    imax=40;                     % Maximum number of iterations for Newton's method 
    Fext = Body.Fext;
    for ii=1:imax    
        tic; 
        
        %[u_bc,deltaf] = Newton_full(Body,Fext);
        [u_bc,deltaf] = Newton_Broyden(ii, Body, Fext); % requires much more steps (~300) and "linear"   

        Body.u(Body.bc) = Body.u(Body.bc)+u_bc;         % Add displacement to previous one
        Body.q(Body.bc) = Body.q(Body.bc)+u_bc;         % change the global positions

        titer=toc;
        titertot=titertot+titer;

        if printStatus(deltaf, u_bc, Re, i, ii, imax, steps, titertot)
            break;  
        end                
    end           

    Body = SaveResults(Body,i,"last");
end

% POST PROCESSING ###############################################
visDeformed = true;
visInitial = true;
PostProcessing(Body,visDeformed,visInitial) 

% Volume change check
faces=Body.BodyFaces;
vertices_before = feval(Body.SurfacefunctionName, Body, Body.q0);
vertices_after  = feval(Body.SurfacefunctionName, Body, Body.q);

V_after = VolumeViaFaces(vertices_after, faces);
V_before = VolumeViaFaces(vertices_before, faces);

fprintf('Volume before: %10.12f; Volume after: %10.12f; Relative change: %10.12f \n', V_before, V_after, (V_after-V_before)/V_before)

CleanTemp(Body, true)

