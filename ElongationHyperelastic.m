clc,clear,close all;
format long

addpath("MainFunctions","MeshFunctions",'InnerForceFunctions',"Postprocessing");

Body.Name = "Body";
CaseName =  string(mfilename);
% ########### Problem data ################################################
Body = DefineElement(Body,"Beam","ANCF",3333,"None");  % 1 - BodyName, 2 - type (beam, plate, etc.), 3 - element name, 4 - modification name (None, EDG, etc.)  
                                                       % ANCF Beam: 3243, 3333, 3343, 3353, 3363, 34X3 (34103)    
[Body,Force,Boundary] = CaseProblemSet(Body,mfilename,"Poigen");  % Itegration Scheme: Poigen, Standard
% ########## Create FE Model ##############################################
ElementNumber = 1;
Body = CreateFEM(Body,ElementNumber);
% % ########## Calculation adjustments ######################################
Body.FiniteDiference= "AceGen"; % Calculation of FD: Matlab, AceGen
Body.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body.DeformationType = "Finite"; % Deformation type: Finite, Small
Body = AddTensors(Body);
% % ########## Visualization of initial situation ##########################
Results = [];  
visualization(Body,Body.q0,'cyan',false); % initial situation
% %####################### Solving ######################################## 
steps = 15;  % sub-loading steps
titertot=0; 

SolutionRegType = "off";  % Regularization type: off, penaltyK, penaltyKf, Tikhonov

Body = CreateBC(Body, Force, Boundary); % Application of Boundary conditions

CreateMex;

%START NEWTON'S METHOD
for i=1:steps
    % Update forces

    Body = SubLoading(Body, i, steps, "linear"); 

    Re=10^(-5);                   % Stopping criterion for residual
    imax=20;                      % Maximum number of iterations for Newton's method 
    Fext = Body.Fext;
    for ii=1:imax    
        tic; 

        [K,Fe] = InnerForce_mex(Body);

        K_bc = K(Body.bc,Body.bc);            % Eliminate linear constraints from stiffness matrix
        ff =  Fe - Fext;

        ff_bc=ff(Body.bc);               % Eliminate linear constraints from force vector
        deltaf=ff_bc/norm(Fext(Body.bc));% Compute residua

        u_bc = Regularization(K_bc,ff_bc,SolutionRegType);  


        if printStatus(deltaf, u_bc, Re, i, ii, imax, steps, titertot)
            break;  
        end                


        Body.u(Body.bc) = Body.u(Body.bc)+u_bc;         % Add displacement to previous one
        Body.q(Body.bc) = Body.q(Body.bc)+u_bc;         % change the global positions

        titer=toc;
        titertot=titertot+titer;   

    end           

    %Pick nodal displacements from result vector
    xlocName = 'xloc' + Body.ElementType;
    uf = Body.u(Body.fextInd); 

    Results = [Results; Body.ElementNumber Body.TotalDofs uf'];
end
% POST PROCESSING ###############################################
visDeformed = false;
visInitial = false;
PostProcessing(Body,Results,visDeformed,visInitial) 
CleanTemp(Body, true)