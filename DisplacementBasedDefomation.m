clc, clear, close all;

addpath("MainFunctions", "MeshFunctions", "InnerForceFunctions","Postprocessing")
Body.Name = "Body";

% ########### Problem data ################################################
% ANCF Beam: 3243, 3333, 3343, 3353, 3363, 34X3 (34103)
Body = DefineElement(Body,"Beam","ANCF",3333,"None");  
Body = Geometry(Body,"ten_Sol_3","Poigen", "Gaus");
Body = Materials(Body,"GOH","optimized_SEE"); 

% ########## Create FE Models #############################################
ElementNumber = 1;
Body = CreateFEM(Body,ElementNumber);

% ########## Calculation adjustments ######################################
Body.FiniteDiference= "AceGen"; % Calculation of FD: Matlab, AceGen
Body.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body.DeformationType = "Finite"; % Deformation type: Finite, Small
Body = AddTensors(Body);

% ########## Boundary Conditions ##########################################
Displacement.Maginutude.X = 0.035;  % Axial displacement appleed to body1 end
Displacement.Maginutude.Y =  0;  
Displacement.Maginutude.Z = 0;  

Displacement.Position.X = Body.Length.X;  % Elongation
Displacement.Position.Y = 0;  
Displacement.Position.Z = 0; 

% Fixed boundaries (applied locally, shift and curvature are accounted automaticaly)
Boundary.Position.X = 0;  
Boundary.Position.Y = 0;
Boundary.Position.Z = 0;
Boundary.Type = "reduced"; % there are s1everal types: full, reduced, positions, none

Body = CreateBCDisplacement(Body, Displacement, Boundary); % Application of Dirichlet boundary conditions
% ####################### Solving ######################################## 
titertot=0;  
Re=1e-4;                   % Stopping criterion for residual
imax=30;                   % Maximum number of iterations for Newton's method 
SolutionRegType = "off";  
Results = [];
bcInd = Body.bcInd; % 0Transformation of body DOF indentifiers to global system

style = "linear";
steps = 10;

%START NEWTON'S METHOD   
for stepnumber=1:steps
    % Prestep
    Body = SubLoadingDispl(Body, stepnumber, steps, style); 
    [K,ff] = InnerForce(Body);% inner forces 

    %% Transformations to accommodate for prescribed displacements        
    ff_bc=ff+K*Body.applied_disp; %modify the whole right-hand side of the system to take into account Dirichlet boundary conditions
    ff_bc(bcInd)=Body.applied_disp(bcInd); %put prescribed values of non-uniform Dirichlet boundary condition
    
    %put zeros in rows and columns corresponding to specified Dirihlet BC
    K_bc = K;
    K_bc(bcInd,:)=0;
    K_bc(:,bcInd)=0;        
    K_bc(sub2ind(size(K),bcInd,bcInd))=-1; %put unity on the main diagonal in the positions corresponding to Dirihlet BC                                       

    % Calculations      
    u_bc = Regularization(K_bc,ff_bc,SolutionRegType); 

    Body.u = Body.u + u_bc(1:Body.TotalDofs);             
    Body.q = Body.q + u_bc(1:Body.TotalDofs);

    titer=toc;
    titertot=titertot+titer;

    for ii=1:imax
        tic;
        % inner forces
        [K,ff] = InnerForce(Body);
        ff_bc = ff(Body.bc);        
        K_bc = K(Body.bc,Body.bc);
                                                    
        % Calculations
        u_bc = Regularization(K_bc,ff_bc,SolutionRegType); 
        Body.u(Body.bc) = Body.u(Body.bc)+u_bc;         % Add displacement to previous one
        Body.q(Body.bc) = Body.q(Body.bc)+u_bc;         % change the global positions
        
        titer=toc;
        titertot=titertot+titer;
    
        if printStatus(ff_bc, u_bc, Re, stepnumber, ii, imax, steps, titertot)            
           break;  
        end 
    end
    uf = Body.u(Body.bcInd);
    Results = [Results; Body.ElementNumber Body.TotalDofs uf(end-2:end)'];    
end 
% POST PROCESSING ###############################################
visualization(Body,Body.q,'cyan',false);
visualization(Body,Body.q0,'blue',false);
PostProcessing(Body,Results,false,false) 
CleanTemp(Body, true)
