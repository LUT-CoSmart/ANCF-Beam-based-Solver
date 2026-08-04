clc,clear,close all;
format long
addpath(genpath(pwd));
Body.Name = "Body";

% ########### Problem data ################################################
Body = DefineElement(Body,"Beam","ANCF",3363,"None");  % 1 - BodyName, 2 - type (beam, plate, etc.), 3 - element name, 4 - modification name (None, EDG, etc.)  
Body = Geometry(Body,'Rectangular',"Standard", "Gauss");  % Cross Sections: Rectangular, Oval, C, Tendon                                                       
Body = Materials(Body,'KS'); % Material models: Gas.-Ogd.-Hol. (GOH), Neo-Hookean (Neo), 2- and 5- constant Mooney-Rivlin (Mooney2, Mooney5),  Kirhhoff-Saint-Venant (KS). 

% Rotation (in degrees)
Body.Rotation.X = 0;
Body.Rotation.Y = 0;
Body.Rotation.Z = 0;

% ########## Boundary Conditions ##########################################
% Force 
Force.Maginutude.X = 0;  % Elongation
Force.Maginutude.Y = 5e8;  
Force.Maginutude.Z = 0;  

Force.Position.X = Body.Length.X / 2;  
Force.Position.Y = 0;  
Force.Position.Z = 0; 

% Boundaries (applied locally, shift and curvature are accounted automaticaly)
Boundary.Position.X = [0 Body.Length.X];  
Boundary.Position.Y = [0 0];
Boundary.Position.Z = [0 0];
Boundary.Type = "full"; % there are several types: full, reduced, positions, none

% ########## Create FE Model ##############################################
ElementNumber = 20;
Body = CreateFEM(Body,ElementNumber,"rectagulars");

% ########## Calculation adjustments ######################################
Body.FiniteDiference= "AceGen"; % Calculation of FD: Matlab, AceGen, Matlab_automatic
Body.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body.DeformationType = "Finite"; % Deformation type: Finite, Small
Body = AddTensors(Body);

% %####################### Solving ######################################## 
steps = 20;  % sub-loading steps
titertot=0;  
Body = CreateBC(Body, Force, Boundary); % Application of Boundary conditions

%START NEWTON'S METHOD   
for i=1:steps
    
    Body = SubLoading(Body, i, steps, "linear"); 
    
    Re=10^(-6);                   % Stopping criterion for residual
    imax=800;                     % Maximum number of iterations for Newton's method 
                                  % it is taken large for Krylov-based (CG) algorithm  
    Fext = Body.Fext;
    for ii=1:imax    
        
        tic;  
        [u_bc,deltaf] = Newton_Broyden(ii, Body, Fext); 
        Body.u(Body.bc) = Body.u(Body.bc)+u_bc;         % Add displacement to previous one
        Body.q(Body.bc) = Body.q(Body.bc)+u_bc;         % change the global positions
        
        titer=toc;
        titertot=titertot+titer;   

        if printStatus(deltaf, u_bc, Re, i, ii, imax, steps, titertot)
            break;  
        end 
       
    end           

    Body = SaveResults(Body,i, "all"); % options: "all", "last", each by (number) 
end

% POST PROCESSING ###############################################
visDeformed = true;
visInitial = false;
PostProcessing(Body,visDeformed,visInitial) 
CleanTemp(Body, true)