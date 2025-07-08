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
Body = Materials(Body,"GOH"); % Material models: GOH (GOH, Amir), Neo-Hookean (Neo), 2- and 5- constant Mooney-Rivlin (Mooney2, Mooney5),  Kirhhoff-Saint-Venant (KS).
% Itegration Scheme: Poigen, Standard
Body = Geometry(Body,"Tendon","Standard");  % Cross Sections: Rectangular, Oval, C, Tendon, Middle_cross_section1_1
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
Body.Twist.angle = 0; % in degrees
Body.Twist.ro = 0;
% ########## Create FE Model ##############################################
ElementNumber = 1;
Body = CreateFEM(Body,ElementNumber);
% ########## Calculation adjustments ######################################
Body.FiniteDiference= "AceGen"; % Calculation of FD: Matlab, AceGen
Body.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body.DeformationType = "Finite"; % Deformation type: Finite, Small
Body = AddTensors(Body);
% ########## Boundary Conditions ##########################################
% Force 
k =1 ;
Force.Maginutude.X = 10*k;  % Elongation
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
Boundary.Type = "reduced"; % there are several types: full, reduced, positions, none

Body = CreateBC(Body, Force, Boundary); % Application of Boundary conditions
% ########## Visualization of initial situation ###########################
Results = [];  
visualization(Body,Body.q0,'red',false);
% % %####################### Solving ######################################## 
steps = 20;  % sub-loading steps
titertot=0;  
Re=10^(-4);           % Stopping criterion for residual
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
    xlocName = 'xloc' + Body.ElementType;
    uf = Body.u(Body.fextInd); 

    Results = [Results; Body.ElementNumber Body.TotalDofs uf'];
end
% POST PROCESSING ###############################################
visDeformed = true;
visInitial = false;
PostProcessing(Body,Results,visDeformed,visInitial) 
CleanTemp(Body, true)

disp = uf(1);
fprintf('Displacement:  %.20f, k = %d\n', disp, k);

% 
% L = Body.Length.Ln;
% H = Body.Length.Y;
% W = Body.Length.Z;
% r1 = Shape_(L,H,W,1,-1,-1)*Body.q;
% r2 = Shape_(L,H,W,-1,-1,-1)*Body.q;
% 
% XX = norm(r1-r2);
% 
% 
% r1 = Shape_(L,H,W,1, 1,-1)*Body.q;
% r2 = Shape_(L,H,W,1, -1,-1)*Body.q;
% YY = norm(r1-r2);
% 
% 
% r1 = Shape_(L,H,W,1, 1, 1)*Body.q;
% r2 = Shape_(L,H,W,1, 1,-1)*Body.q;
% ZZ = norm(r1-r2);
% 
% % Body.Volume
% % 
% % VV = XX*ZZ*YY