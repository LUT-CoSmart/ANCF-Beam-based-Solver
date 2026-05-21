clc,clear,close all;
format long
addpath(genpath(pwd));
Body.Name = "Body";
CaseName =  string(mfilename);
CaseSubtype = "Large"; % there are two options: Large & Small

% ########### Problem data ################################################
Body = DefineElement(Body,"Beam","ANCF",3363,"None");  % 1 - BodyName, 2 - type (beam, plate, etc.), 3 - element name, 4 - modification name (None, EDG, etc.)  
                                                       % ANCF Beam: 3243, 3333, 3343, 3353, 3363, 34X3 (34103)    
[Body,Force,Boundary] = CaseProblemSet(Body,mfilename + CaseSubtype,"Standard");  % Itegration Scheme: Poigen, Standard

% ########## Create FE Model ##############################################
ElementNumber = 16;
Body = CreateFEM(Body,ElementNumber);

% ########## Calculation adjustments ######################################
Body.FiniteDiference= "AceGen"; % Calculation of FD: Matlab, AceGen
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
    
    Re=10^(-4);                   % Stopping criterion for residual
    imax=800;                     % Maximum number of iterations for Newton's method 
                                  % it is taken large for Krylov-based (CG) algorithm  
    Fext = Body.Fext;
    for ii=1:imax    
        
        tic; 
        % [u_bc,deltaf] = Newton_full(Body,Fext);  
        [u_bc,deltaf] = Newton_Broyden(ii, Body, Fext); 
        % [u_bc,deltaf] = Newton_BFGS(ii, Body, Fext)
        % [u_bc,deltaf] = Newton_Krylov(ii, Body, Fext, Re, "JF"); % options: CG - Conjugate Gradient, JF - Jacobian Free  
        
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
visInitial = true;
PostProcessing(Body,visDeformed,visInitial) 
% CleanTemp(Body, true)