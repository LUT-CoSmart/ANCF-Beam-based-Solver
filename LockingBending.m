clc,clear,close all;
format long
addpath("MainFunctions");
addpath("MeshFunctions");
addpath("Postprocessing");
Body.Name = "Body";
CaseName =  string(mfilename);
CaseSubtype = "Large"; % there are two options: Large & Small
% ########### Problem data ################################################
Body = DefineElement(Body,"Beam","ANCF",3333,"None");  % 1 - BodyName, 2 - type (beam, plate, etc.), 3 - element name, 4 - modification name (None, EDG, etc.)  
                                                       % ANCF Beam: 3243, 3333, 3343, 3353, 3363, 34X3 (34103)    
[Body,Force,Boundary] = CaseProblemSet(Body,mfilename + CaseSubtype,"Standard");  % Itegration Scheme: Poigen, Standard
% ########## Create FE Model ##############################################
ElementNumber = 4;
Body = CreateFEM(Body,ElementNumber);
% ########## Calculation adjustments ######################################
Body.FiniteDiference= "AceGen"; % Calculation of FD: Matlab, AceGen
Body.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body.DeformationType = "Finite"; % Deformation type: Finite, Small
Body = AddTensors(Body);
% ########## Visualization of initial situation ###########################
Results = [];  
% visualization(Body,Body.q0,'cyan',false); % initial situation
% %####################### Solving ######################################## 
steps = 2;  % sub-loading steps
titertot=0;  
%START NEWTON'S METHOD   
for i=1:steps
    
    % Update forces, supported loading types: linear, exponential, quadratic, cubic;
    Subforce = SubLoading(Force, i, steps, "linear"); 

    % Application of Boundary conditions
    Body = CreateBC(Body, Subforce, Boundary);
    
    Re=10^(-4);                   % Stopping criterion for residual
    imax=20;                      % Maximum number of iterations for Newton's method 
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
            fprintf('Convergence: %10.4f, Displacements norm: %10.5f\n', norm(abs(deltaf)), norm(u_bc));
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