clc,clear,close all;
format long
addpath("MainFunctions", "MeshFunctions", "InnerForceFunctions","Postprocessing")
addpath(genpath("Contact"))

Body1.Name = "Body1";
Body2.Name = "Body2";
Body3.Name = "Body3";
% ########### Problem data ################################################
% ANCF Beam: 3243, 3333, 3343, 3353, 3363, 34X3 (34103)
Body1 = DefineElement(Body1,"Beam","ANCF",3363,"None");  
Body2 = DefineElement(Body2,"Beam","ANCF",3363,"None"); 
Body3 = DefineElement(Body3,"Beam","ANCF",3363,"None"); 
% Material models: GOH (GOH), Neo-Hookean (Neo), 2- and 5- constant Mooney-Rivlin (Mooney2, Mooney5),  Kirhhoff-Saint-Venant (KS).
Body1 = Materials(Body1,'KS'); 
Body2 = Materials(Body2,'KS'); 
Body3 = Materials(Body3,'KS'); 
% Geometry
Body1 = Geometry(Body1,"Rectangular","Standard");  % Cross Sections: Rectangular, Oval, C, Tendon
Body2 = Geometry(Body2,"Rectangular","Standard");  % Itegration Scheme: Poigen, Standard
Body3 = Geometry(Body3,"Rectangular","Standard");  % Itegration Scheme: Poigen, Standard
% ########### Set Bodies positions ########################################
% Shift of Body1
Body1.Shift.X =  0;
Body1.Shift.Y =  Body2.Length.Y;
Body1.Shift.Z =  0;

Body3.Shift.X =  0;
Body3.Shift.Y = -Body2.Length.Y;
Body3.Shift.Z =  0;

% ########## Create FE Models #############################################

ElementNumber1 = 2;
Body1 = CreateFEM(Body1,ElementNumber1);
ElementNumber2 = 2;
Body2 = CreateFEM(Body2,ElementNumber2);
ElementNumber3 = 2;
Body3 = CreateFEM(Body3,ElementNumber2);

% ########## Calculation adjustments ######################################
Body1.FiniteDiference= "Matlab"; % Calculation of FD: Matlab, AceGen
Body1.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body1.DeformationType = "Finite"; % Deformation type: Finite, Small
Body1 = AddTensors(Body1);

Body2.FiniteDiference= "Matlab"; % Calculation of FD: Matlab, AceGen
Body2.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body2.DeformationType = "Finite"; % Deformation type: Finite, Small
Body2 = AddTensors(Body2);

Body3.FiniteDiference= "Matlab"; % Calculation of FD: Matlab, AceGen
Body3.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body3.DeformationType = "Finite"; % Deformation type: Finite, Small
Body3 = AddTensors(Body3);

% ########## Boundary Conditions ##########################################
% Body1 
% Force (applied locally, shift and curvature are accounted automaticaly)
Force1.Maginutude.Y = -62.5*10^(5);  
Force1.Position.X = Body1.Length.X;  % Elongation

% Boundaries (applied locally, shift and curvature are accounted automaticaly)
Boundary1.Position = [];
Boundary1.Type = "full"; % there are several types: full, reduced, positions, none

% Body2
Force2.Maginutude = [];
Force2.Position.X = Body2.Length.X;  % Elongation

% Boundaries
Boundary2.Position = [];
Boundary2.Type = "full"; % there are several types: full, reduced, positions, none


% Body2
Force3.Maginutude = [];
Force3.Position.X = Body3.Length.X;  % Elongation

% Boundaries
Boundary3.Position = [];
Boundary3.Type = "full"; % there are several types: full, reduced, positions, none

% ########## Contact characteristics ######################################
ContactType = "NitscheLin"; % Options: "None", "Penalty", "NitscheLin"...
ContactVariable = 1e9;
Body1.ContactRole = "slave"; % Options: "master", "slave"
Body2.ContactRole = "master";
Body3.ContactRole = "slave";
% ########## Visualization of initial situation ###########################
% figure;
% hold on
% axis equal 
% xlabel('\it{X}','FontName','Times New Roman','FontSize',[20])
% ylabel('\it{Y}','FontName','Times New Roman','FontSize',[20]),
% zlabel('Z [m]','FontName','Times New Roman','FontSize',[20]);
% visualization(Body1,Body1.q0,'cyan',true);
% visualization(Body2,Body2.q0,'red',true);

% %####################### Solving ######################################## 
steps = 10;  % sub-loading steps
titertot=0;  
Re=10^(-3);                   % Stopping criterion for residual
imax=15;                      % Maximum number of iterations for Newton's method 
SolutionRegType = "off";  % Regularization type: off, penaltyK, penaltyKf, Tikhonov
ContactRegType = "off";
Results1 = [];
Results2 = [];
Results3 = [];

Body1 = CreateBC(Body1, Force1, Boundary1); % Application of Boundary conditions
Body2 = CreateBC(Body2, Force2, Boundary2); % Application of Boundary conditions
Body3 = CreateBC(Body3, Force3, Boundary3); % Application of Boundary conditions

%START NEWTON'S METHOD   
for i=1:steps

    Body1 = SubLoading(Body1, i, steps, "linear"); 
    Body2 = SubLoading(Body2, i, steps, "cubic"); 
    Body3 = SubLoading(Body3, i, steps, "cubic"); 

    Fext1 = Body1.Fext;
    Fext2 = Body2.Fext;
    Fext3 = Body3.Fext;


    bc = [Body1.bc Body2.bc Body3.bc];

    for ii=1:imax
        tic;

        % Contact forces
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Kc1,Fc1,Gap1] = Contact(Body1,Body2,ContactType,ContactVariable,ContactRegType);
        Fc1_extend = [Fc1; zeros(Body3.TotalDofs,1)];
        Kc1_extend = [Kc1 zeros(Body1.TotalDofs+Body2.TotalDofs, Body3.TotalDofs);
                      zeros(Body3.TotalDofs,Body1.TotalDofs+Body2.TotalDofs + Body3.TotalDofs)];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Kc2,Fc2,Gap2] = Contact(Body2,Body3,ContactType,ContactVariable,ContactRegType);
        Fc2_extend = [zeros(Body1.TotalDofs,1); Fc2];

        Kc2_extend = [zeros(Body1.TotalDofs, Body1.TotalDofs + Body2.TotalDofs+Body3.TotalDofs);
               zeros(Body2.TotalDofs+Body3.TotalDofs, Body1.TotalDofs) Kc2];    
        
        Fc = Fc1_extend + Fc2_extend;     
        Kc = Kc1_extend + Kc2_extend;
        
        Gap = Gap1 + Gap2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % inner forces
        [Ke1,Fe1] = InnerForce(Body1);
        [Ke2,Fe2] = InnerForce(Body2);
        [Ke3,Fe3] = InnerForce(Body3);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assemblance
        Fe = [Fe1; Fe2 ; Fe3];
        Fext = [Fext1; Fext2 ; Fext3];
        Ke = [Ke1 zeros(Body1.TotalDofs,Body2.TotalDofs+ Body3.TotalDofs);
              zeros(Body2.TotalDofs,Body1.TotalDofs) Ke2 zeros(Body2.TotalDofs,Body3.TotalDofs);
              zeros(Body3.TotalDofs,Body1.TotalDofs+ Body2.TotalDofs) Ke3];

        K = Kc + Ke;
        ff =  Fe - Fext + Fc;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculations
        K_bc = K(bc,bc); 
        ff_bc = ff(bc);
        deltaf=ff_bc/norm(Fext(bc)); 
        u_bc = Regularization(K_bc,ff_bc,SolutionRegType); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Separation
        Body1.u(Body1.bc) = Body1.u(Body1.bc) + u_bc(1:Body1.ndof);
        Body2.u(Body2.bc) = Body2.u(Body2.bc) + u_bc(Body1.ndof + 1:Body1.ndof + Body2.ndof);
        Body3.u(Body3.bc) = Body3.u(Body3.bc) + u_bc(Body1.ndof + Body2.ndof + 1:end);
        
        Body1.q(Body1.bc) = Body1.q(Body1.bc) + u_bc(1:Body1.ndof);
        Body2.q(Body2.bc) = Body2.q(Body2.bc) + u_bc(Body1.ndof + 1:Body1.ndof + Body2.ndof);
        Body3.q(Body3.bc) = Body3.q(Body3.bc) + u_bc(Body1.ndof + Body2.ndof + 1:end);

        titer=toc;
        titertot=titertot+titer;

        if printStatus(deltaf, u_bc, Re, i, ii, imax, steps, titertot, Gap)
            break;  
        end 

    end
   


    %Pick nodal displacements from result vector
    xlocName1 = 'xloc' + Body1.ElementType;
    uf1 = Body1.u(Body1.fextInd);
    Results1 = [Results1; Body1.ElementNumber Body1.TotalDofs uf1'];

    xlocName2 = 'xloc' + Body2.ElementType;
    uf2 = Body2.u(Body2.fextInd); 
    Results2 = [Results2; Body2.ElementNumber Body2.TotalDofs uf2'];

    xlocName3 = 'xloc' + Body3.ElementType;
    uf3 = Body3.u(Body3.fextInd); 
    Results3 = [Results3; Body3.ElementNumber Body3.TotalDofs uf3'];
end    

% POST PROCESSING ###############################################
hold on
xlabel('\it{X}','FontName','Times New Roman','FontSize',[20])
ylabel('\it{Y}','FontName','Times New Roman','FontSize',[20]),
zlabel('Z [m]','FontName','Times New Roman','FontSize',[20]);
visualization(Body1,Body1.q,'cyan',true);
visualization(Body2,Body2.q,'none',true);
visualization(Body3,Body3.q,'blue',true);

PostProcessing(Body1,Results1,false,false) 
PostProcessing(Body2,Results2,false,false) 
PostProcessing(Body3,Results3,false,false)

CleanTemp(Body1, true)
CleanTemp(Body2, true)
CleanTemp(Body3, true)