clc,clear,close all;
format long
addpath("MainFunctions");
addpath("Postprocessing");
addpath("MeshFunctions")
CaseName =  string(mfilename);
Body.Name = "Body";
%% There are two options: Large & Small
CaseSubtype = "Large"; 
% ########### Problem data ################################################
Body = DefineElement(Body,"Beam","ANCF",3333,"None");  % 1 - BodyName, 2 - type (beam, plate, etc.), 3 - element name, 4 - modification name (None, EDG, etc.)  
                                                       % ANCF Beam: 3243, 3333, 3343, 3353, 3363, 34X3 (34103)    
[Body,Force,Boundary] = CaseProblemSet(Body,mfilename + CaseSubtype,"Standard");  % Itegration Scheme: Poigen, Standard
% ########## Create FE Model ##############################################
ElementNumber = 1;
Body = CreateFEM(Body,ElementNumber);
% ########## Calculation adjustments ######################################
Body.FiniteDiference= "AceGen"; % Calculation of FD: Matlab, AceGen
Body.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body.DeformationType = "Finite"; % Deformation type: Finite, Small
Body = AddTensors(Body);
% ########## Visualization of initial situation ###########################
Results = [];  
visualization(Body,Body.q0,'cyan',false); % initial situation
% %####################### Solving ########################################
angleSet = 0:5:90; 
steps = 1;  % sub-loading steps
SolutionRegType = "off";  % Regularization type: off, penaltyK, penaltyKf, Tikhonov

for angle = angleSet
    fprintf("Considered angle %d\n", angle);
    titertot=0; 

    Force.Maginutude.Y = -Force.Maginutude.Total * cosd(angle);
    Force.Maginutude.Z =  Force.Maginutude.Total * sind(angle);
     
    Body = CreateBC(Body, Force, Boundary); % Application of Boundary conditions

    %START NEWTON'S METHOD   
    for i=1:steps
        
        Body = SubLoading(Body, i, steps, "linear");    
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

    end
    Results = [Results; Body.ElementNumber Body.TotalDofs uf'];
end    
% POST PROCESSING ###############################################
visDeformed = false;
visInitial = false;
PostProcessing(Body,Results,visDeformed,visInitial) 


subplot(2,1,1);
plot(angleSet, -Results(:,4), 'b--');
title('Displacement along Y');
xlabel('Angle');
ylabel('Y');
grid on;

subplot(2,1,2);
plot(angleSet, Results(:,5), 'r--');
title('Displacement along Z');
xlabel('Angle');
ylabel('Z');
grid on;

CleanTemp(Body, true)