clc,clear,close all;
format long

addpath("MainFunctions");
addpath("MeshFunctions")
addpath("Postprocessing");
addpath('InnerForceFunctions');

Body.Name = "Body";
% ########### Problem data ################################################
Body = DefineElement(Body,"Beam","ANCF",3333,"None");  % 1 - BodyName, 2 - type (beam, plate, etc.), 3 - element name, 4 - modification name (None, EDG, etc.)  
                                                       % ANCF Beam: 3243, 3333, 3343, 3353, 3363, 34X3 (34103)    
Body = Materials(Body,"GOH", "Amir"); % Material models: GOH (GOH, Amir), Neo-Hookean (Neo), 2- and 5- constant Mooney-Rivlin (Mooney2, Mooney5),  Kirhhoff-Saint-Venant (KS).
                                                        % Integration Scheme: Poigen, Standard
Body = Geometry(Body,"ten_Sol_3","Poigen", "Lobatto");  % Cross Sections: Rectangular, Oval, C, Tendon, etc.
                                                        % Integration points of generating line: Gauss, Lobatto
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
Body.FiniteDiference= "Matlab_automatic"; % Calculation of FD: Matlab, Matlab_automatic, AceGen
Body.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body.DeformationType = "Finite"; % Deformation type: Finite, Small
Body = AddTensors(Body);
% ########## Boundary Conditions ##########################################
% Force 
steps = 40;  % sub-loading steps
Force.Maginutude.X = 3e+04;  % Elongation

% steps = 50;  % sub-loading steps
% Force.Maginutude.X = 2.889627761315940*1e5;
% Force.Maginutude.Y = 0;  
% Force.Maginutude.Z = 0;  

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
Results = [];  
visualization(Body,Body.q0,'red',false);
% % %####################### Solving ######################################## 
%steps = 40;  % sub-loading steps
titertot=0;  
Re=1*10^(-4);           % Stopping criterion for residual
imax=40;              % Maximum number of iterations for Newton's method 

% origFolder = pwd;
% cd InnerForceFunctions;
% codegen InnerForce -args {Body} -config:mex
% cd(origFolder);

SolutionRegType = "off";      % Regularization type: off, penaltyK, penaltyKf, Tikhonov
%START NEWTON'S METHOD   

for i=1:steps

    % Update forces, supported loading types: linear, exponential, quadratic, cubic;
    Body = SubLoading(Body, i, steps, "cubic"); 
    Fext = Body.Fext;    
                    
    for ii=1:imax    
        tic; 

        [K,Fe] = InnerForce(Body);
                       
        K_bc = K(Body.bc,Body.bc);            % Eliminate linear constraints from stiffness matrix
        ff =  Fe - Fext;

        ff_bc=ff(Body.bc);               % Eliminate linear constraints from force vector
        deltaf=ff_bc/norm(Fext(Body.bc));% Compute residual
        
        u_bc = Regularization(K_bc,ff_bc,SolutionRegType);  % Compute displacements   

        Body.u(Body.bc) = Body.u(Body.bc)+u_bc;         % Add displacement to previous one
        Body.q(Body.bc) = Body.q(Body.bc)+u_bc;         % change the global positions
        titer=toc;
        titertot=titertot+titer;   


        if printStatus(deltaf, u_bc, Re, i, ii, imax, steps, titertot)
            break;  
        end 
      

        
    end           

    %Pick nodal displacements from result vector
    uf = Body.u(Body.fextInd); 
    Results = [Results; Body.ElementNumber Body.TotalDofs uf'];
    
end
% POST PROCESSING ###############################################
visDeformed = true;
visInitial = true;
PostProcessing(Body,Results,visDeformed,visInitial) 
% CleanTemp(Body, true)

% volume check
faces=Body.BodyFaces;
vertices_before = feval(Body.SurfacefunctionName, Body, Body.q0);
vertices_after  = feval(Body.SurfacefunctionName, Body, Body.q);

V_after = VolumeViaFaces(vertices_after, faces);
V_before = VolumeViaFaces(vertices_before, faces);

fprintf('Volume before: %10.12f; Volume after: %10.12f; Relative change: %10.12f \n', V_before, V_after, (V_after-V_before)/V_before)
