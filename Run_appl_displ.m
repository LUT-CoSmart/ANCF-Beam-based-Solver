clc,clear,close all;
format long
addpath("MainFunctions")
addpath("MeshFunctions")
addpath("InnerForceFunctions")
addpath(genpath("Contact"))
addpath("Postprocessing");
Body1.Name = "Body1";
Body2.Name = "Body2";
Body3.Name = "Body3";
% ########### Problem data ################################################
% ANCF Beam: 3243, 3333, 3343, 3353, 3363, 34X3 (34103)
Body1 = DefineElement(Body1,"Beam","ANCF",3333,"None");  
Body2 = DefineElement(Body2,"Beam","ANCF",3333,"None");  
Body3 = DefineElement(Body3,"Beam","ANCF",3333,"None"); 
% Material models: GOH (GOH), Neo-Hookean (Neo), 2- and 5- constant Mooney-Rivlin (Mooney2, Mooney5),  Kirhhoff-Saint-Venant (KS).
Body1 = Materials(Body1,"Neo","Alex"); 
Body2 = Materials(Body2,"Neo","Alex"); 
Body3 = Materials(Body3,'Neo',"Alex");
% Geometry
Body1 = Geometry(Body1,"MG_subj2_middle","Poigen");  % Cross Sections: Rectangular, Oval, C, Tendon
Body2 = Geometry(Body2,"Sol_subj2_middle","Poigen");  % Itegration Scheme: Poigen, Standard
Body3 = Geometry(Body3,"LG_subj2_middle","Poigen");  % Itegration Scheme: Poigen, Standard
% ########### Set Bodies positions ########################################
angle = 2;
% Tendon twist
Center1 = [Body1.CSCenterY, Body1.CSCenterZ];
Center2 = [Body2.CSCenterY, Body2.CSCenterZ];
Center3 = [Body3.CSCenterY, Body3.CSCenterZ];

RelCenter1 = Center1 - Center2;
RelCenter2 = [0, 0];
RelCenter3 = Center3 - Center2;

Body1.Twist.angle = angle;
Body1.Twist.initial_rot = atan2d(RelCenter1(2),RelCenter1(1));
Body1.Twist.ro = norm(Center2 - Center1);

Body2.Twist.angle = angle;
Body2.Twist.initial_rot = atan2d(RelCenter2(2),RelCenter2(1));
Body2.Twist.ro = 0;

Body3.Twist.angle = angle;
Body3.Twist.initial_rot = atan2d(RelCenter3(2),RelCenter3(1));
Body3.Twist.ro = norm(Center2 - Center3);

% Rotation (in degrees)
Body1.Rotation.X = 0;
Body1.Rotation.Y = 0;
Body1.Rotation.Z = 0;

Body2.Rotation.X = 0;
Body2.Rotation.Y = 0;
Body2.Rotation.Z = 0;

Body3.Rotation.X = 0;
Body3.Rotation.Y = 0;
Body3.Rotation.Z = 0;

% ########## Create FE Models #############################################

ElementNumber1 = 4;
Body1 = CreateFEM(Body1,ElementNumber1);
ElementNumber2 = 4;
Body2 = CreateFEM(Body2,ElementNumber2);
ElementNumber3 = 4;
Body3 = CreateFEM(Body3,ElementNumber3);

% ########## Calculation adjustments ######################################
Body1.FiniteDiference= "AceGen"; % Calculation of FD: Matlab, AceGen
Body1.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body1.DeformationType = "Finite"; % Deformation type: Finite, Small
Body1 = AddTensors(Body1);

Body2.FiniteDiference= "AceGen"; % Calculation of FD: Matlab, AceGen
Body2.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body2.DeformationType = "Finite"; % Deformation type: Finite, Small
Body2 = AddTensors(Body2);

Body3.FiniteDiference= "AceGen"; % Calculation of FD: Matlab, AceGen
Body3.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body3.DeformationType = "Finite"; % Deformation type: Finite, Small
Body3 = AddTensors(Body3);


totalArea = Body1.detF0*sum(Body1.Gint(:,4)) + Body2.detF0*sum(Body2.Gint(:,4)) + Body3.detF0*sum(Body3.Gint(:,4));
frac1 = Body1.detF0*sum(Body1.Gint(:,4)) / totalArea;
frac2 = Body2.detF0*sum(Body2.Gint(:,4)) / totalArea;
frac3 = Body3.detF0*sum(Body3.Gint(:,4)) / totalArea;

% ########## Boundary Conditions ##########################################

Force = 0; % Dirichlet boundary conditions

% Body1 
% Force (applied locally, shift and curvature are accounted automaticaly)
Displacement1.Maginutude.X =  10E-3;  % Axial displacement appleed to body1 end
Displacement1.Maginutude.Y = 0;  
Displacement1.Maginutude.Z = 0;  

Displacement1.Position.X = Body1.Length.X;  % Elongation
Displacement1.Position.Y = 0;  
Displacement1.Position.Z = 0; 

Force1.Maginutude.X =  Force/3;  % Elongation
Force1.Maginutude.Y = 0;  
Force1.Maginutude.Z = 0;  

Force1.Position.X = Body1.Length.X;  % Elongation
Force1.Position.Y = 0;  
Force1.Position.Z = 0; 

% Fixed boundaries (applied locally, shift and curvature are accounted automaticaly)
Boundary1.Position.X = 0;  
Boundary1.Position.Y = 0;
Boundary1.Position.Z = 0;

Boundary1.Type = "reduced"; % there are s1everal types: full, reduced, positions, none

% Body2
%Force2.Maginutude.X = Force*frac2;  % Elongation
Displacement2.Maginutude.X = 2E-3;
Displacement2.Maginutude.Y = 0;  
Displacement2.Maginutude.Z = 0;  

Displacement2.Position.X = Body2.Length.X;  % Elongation
Displacement2.Position.Y = 0;  
Displacement2.Position.Z = 0; 

Force2.Maginutude.X = Force/3;
Force2.Maginutude.Y = 0;  
Force2.Maginutude.Z = 0;  

Force2.Position.X = Body2.Length.X;  % Elongation
Force2.Position.Y = 0;  
Force2.Position.Z = 0; 

% Boundaries
%Boundary2.Position = [];
Boundary2.Position.X = 0;  
Boundary2.Position.Y = 0;
Boundary2.Position.Z = 0;

Boundary2.Type = "reduced"; % there are several types: full, reduced, positions, none

% Body3
%Force3.Maginutude.X = Force*frac3;  % Elongation
Displacement3.Maginutude.X = 5E-3;  % Elongation
Displacement3.Maginutude.Y = 0;  
Displacement3.Maginutude.Z = 0;  
 
Displacement3.Position.X = Body3.Length.X;  % Elongation
Displacement3.Position.Y = 0;  
Displacement3.Position.Z = 0; 

Force3.Maginutude.X = Force/3;  % Elongation
Force3.Maginutude.Y = 0;  
Force3.Maginutude.Z = 0;  
 
Force3.Position.X = Body3.Length.X;  % Elongation
Force3.Position.Y = 0;  
Force3.Position.Z = 0; 

% Fixed Boundaries
%Boundary2.Position = [];
Boundary3.Position.X = 0;  
Boundary3.Position.Y = 0;
Boundary3.Position.Z = 0;

Boundary3.Type = "reduced"; % there are several types: full, reduced, positions, none

% ########## Contact characteristics ######################################
ContactType = "Penalty"; % Options: "None", "Penalty", "NitscheLin"...
ContactVariable = 1e1;
Body1.ContactRole = "master"; % Options: "master", "slave"
Body2.ContactRole = "slave";
Body3.ContactRole = "master";

% %####################### Solving ######################################## 
steps = 20;  % sub-loading steps
titertot=0;  
Re=10^(-4);                   % Stopping criterion for residual
imax=20;                      % Maximum number of iterations for Newton's method 
SolutionRegType = "off";  % Regularization type: off, penaltyK, penaltyKf, Tikhonov
ContactRegType = "off";
Results1 = [];
Results2 = [];
Results3 = [];

%Body1 = CreateBC(Body1, Force1, Boundary1); % Application of Boundary conditions
%Body2 = CreateBC(Body2, Force2, Boundary2); % Application of Boundary conditions
%Body3 = CreateBC(Body3, Force3, Boundary3); % Application of Boundary conditions

Body1 = CreateAllBC(Body1, Force1, Displacement1, Boundary1); % Application of Dirichlet boundary conditions
Body2 = CreateAllBC(Body2, Force2, Displacement2, Boundary2); % Application of Dirichlet boundary conditions
Body3 = CreateAllBC(Body3, Force3, Displacement3, Boundary3); % Application of Dirichlet boundary conditions

% Transformation of body DOF indentifiers to global system

Body2bcIndGlobal = Body2.bcInd + Body1.TotalDofs;
Body3bcIndGlobal = Body3.bcInd + Body1.TotalDofs+ Body2.TotalDofs;

Applied_disp = [Body1.applied_disp Body2.applied_disp Body3.applied_disp];
Applied_zeros = zeros(1,Body1.TotalDofs+Body2.TotalDofs+Body2.TotalDofs);
bcInd = [Body1.bcInd Body2bcIndGlobal Body3bcIndGlobal];

%style = "cubic";
%style = "linear";
%START NEWTON'S METHOD   
%for i=1:steps
    
    % Body1 = SubLoading(Body1, i, steps, style); 
    % Body2 = SubLoading(Body2, i, steps, style); 
    % Body3 = SubLoading(Body3, i, steps, style); 

 
    Fext1 = Body1.Fext;
    Fext2 = Body2.Fext;
    Fext3 = Body3.Fext;

    bc = [Body1.bc Body2.bc Body3.bc];


    % First iteration with non-zero Dirichlet boundary conditions   

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Contact forces
        [Kc1,Fc1,Gap1] = Contact(Body1,Body2,ContactType,ContactVariable,ContactRegType);
        Fc1_extend = [Fc1; zeros(Body3.TotalDofs,1)];
        Kc1_extend = [Kc1 zeros(Body1.TotalDofs+Body2.TotalDofs, Body3.TotalDofs);
                      zeros(Body3.TotalDofs,Body1.TotalDofs+Body2.TotalDofs + Body3.TotalDofs)];
        
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Contact forces
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

        %% Transformations to accommodate for prescribed displacements
        
        ff_bc=ff+K*Applied_disp.'; %modify the whole right-hand side of the system to take into account Dirichlet boundary conditions
        ff_bc(bcInd)=Applied_disp(bcInd).'; %put prescribed values of non-uniform Dirichlet boundary conditions into the correspondin
        K_bc = K;
        %put zeros in rows and columns corresponding to specified Dirihlet BC
        K_bc(bcInd,1:(Body1.TotalDofs+Body2.TotalDofs+Body3.TotalDofs))=0;
        K_bc(1:(Body1.TotalDofs+Body2.TotalDofs+Body3.TotalDofs),bcInd)=0;
        %put unity on the main diagonal in the positions corresponding to Dirihlet BC
        bc_lin = sub2ind(size(K),bcInd,bcInd);   
        K_bc(bc_lin)=-1;                                        


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculations
        %K_bc = K(bc,bc); 
        %ff_bc = ff(bc);

        %deltaf=ff_bc/norm(Fext(bc)); 
        %u_bc = Regularization(K_bc,ff_bc,SolutionRegType); 
        deltaf=ff_bc/norm(Fext); 
        u_bc = Regularization(K_bc,ff_bc,SolutionRegType); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Separation
        % Body1.u(Body1.bc) = Body1.u(Body1.bc) + u_bc(1:Body1.ndof);
        % Body2.u(Body2.bc) = Body2.u(Body2.bc) + u_bc(Body1.ndof + 1:Body1.ndof + Body2.ndof);
        % Body3.u(Body3.bc) = Body3.u(Body3.bc) + u_bc(Body1.ndof + Body2.ndof + 1:end);
        % 
        % Body1.q(Body1.bc) = Body1.q(Body1.bc) + u_bc(1:Body1.ndof);
        % Body2.q(Body2.bc) = Body2.q(Body2.bc) + u_bc(Body1.ndof + 1:Body1.ndof + Body2.ndof);
        % Body3.q(Body3.bc) = Body3.q(Body3.bc) + u_bc(Body1.ndof + Body2.ndof + 1:end);
        Body1.u = Body1.u + u_bc(1:Body1.TotalDofs);
        Body2.u = Body2.u + u_bc(Body1.TotalDofs + 1:Body1.TotalDofs + Body2.TotalDofs);
        Body3.u = Body3.u + u_bc(Body1.TotalDofs + Body2.TotalDofs + 1:end);
        
        Body1.q = Body1.q + u_bc(1:Body1.TotalDofs);
        Body2.q = Body2.q + u_bc(Body1.TotalDofs + 1:Body1.TotalDofs + Body2.TotalDofs);
        Body3.q = Body3.q + u_bc(Body1.TotalDofs + Body2.TotalDofs + 1:end);

        titer=toc;
        titertot=titertot+titer;

    for ii=1:imax
        tic;
        
       % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %  % Contact forces
       %  [Kc1,Fc1,Gap1] = Contact(Body1,Body2,ContactType,ContactVariable,ContactRegType);
       %  Fc1_extend = [Fc1; zeros(Body3.TotalDofs,1)];
       %  Kc1_extend = [Kc1 zeros(Body1.TotalDofs+Body2.TotalDofs, Body3.TotalDofs);
       %                zeros(Body3.TotalDofs,Body1.TotalDofs+Body2.TotalDofs + Body3.TotalDofs)];
       % 
       % 
       %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %  % Contact forces
       %  [Kc2,Fc2,Gap2] = Contact(Body2,Body3,ContactType,ContactVariable,ContactRegType);
       %  Fc2_extend = [zeros(Body1.TotalDofs,1); Fc2];
       % 
       %  Kc2_extend = [zeros(Body1.TotalDofs, Body1.TotalDofs + Body2.TotalDofs+Body3.TotalDofs);
       %         zeros(Body2.TotalDofs+Body3.TotalDofs, Body1.TotalDofs) Kc2];    
       % 
       %  Fc = Fc1_extend + Fc2_extend;     
       %  Kc = Kc1_extend + Kc2_extend;
       %  Gap = Gap1 + Gap2;
       % 
       %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %  % inner forces
       %  [Ke1,Fe1] = InnerForce(Body1);
       %  [Ke2,Fe2] = InnerForce(Body2);
       %  [Ke3,Fe3] = InnerForce(Body3);
       % 
       %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %  % Assemblance
       %  Fe = [Fe1; Fe2 ; Fe3];
       %  Fext = [Fext1; Fext2 ; Fext3];
       %  Ke = [Ke1 zeros(Body1.TotalDofs,Body2.TotalDofs+ Body3.TotalDofs);
       %        zeros(Body2.TotalDofs,Body1.TotalDofs) Ke2 zeros(Body2.TotalDofs,Body3.TotalDofs);
       %        zeros(Body3.TotalDofs,Body1.TotalDofs+ Body2.TotalDofs) Ke3];
       % 
       %  K = Kc + Ke;
       % 
       %  ff =  Fe - Fext + Fc;
       % 
       %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %  % Calculations
       %  K_bc = K(bc,bc); 
       %  ff_bc = ff(bc);
       %  deltaf=ff_bc/norm(Fext(bc)); 
       %  u_bc = Regularization(K_bc,ff_bc,SolutionRegType); 
       % 
       %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %  % Separation
       %  Body1.u(Body1.bc) = Body1.u(Body1.bc) + u_bc(1:Body1.ndof);
       %  Body2.u(Body2.bc) = Body2.u(Body2.bc) + u_bc(Body1.ndof + 1:Body1.ndof + Body2.ndof);
       %  Body3.u(Body3.bc) = Body3.u(Body3.bc) + u_bc(Body1.ndof + Body2.ndof + 1:end);
       % 
       %  Body1.q(Body1.bc) = Body1.q(Body1.bc) + u_bc(1:Body1.ndof);
       %  Body2.q(Body2.bc) = Body2.q(Body2.bc) + u_bc(Body1.ndof + 1:Body1.ndof + Body2.ndof);
       %  Body3.q(Body3.bc) = Body3.q(Body3.bc) + u_bc(Body1.ndof + Body2.ndof + 1:end);
       % 
       %  titer=toc;
       %  titertot=titertot+titer;


           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Contact forces
        [Kc1,Fc1,Gap1] = Contact(Body1,Body2,ContactType,ContactVariable,ContactRegType);
        Fc1_extend = [Fc1; zeros(Body3.TotalDofs,1)];
        Kc1_extend = [Kc1 zeros(Body1.TotalDofs+Body2.TotalDofs, Body3.TotalDofs);
                      zeros(Body3.TotalDofs,Body1.TotalDofs+Body2.TotalDofs + Body3.TotalDofs)];
        
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Contact forces
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

        %% Transformations to accommodate for prescribed displacements
        
        ff_bc=ff; 
        ff_bc(bcInd)=Applied_zeros(bcInd).'; %put prescribed values of non-uniform Dirichlet boundary conditions into the correspondin
        K_bc = K;
        %put zeros in rows and columns corresponding to specified Dirihlet BC
        K_bc(bcInd,1:(Body1.TotalDofs+Body2.TotalDofs+Body3.TotalDofs))=0;
        K_bc(1:(Body1.TotalDofs+Body2.TotalDofs+Body3.TotalDofs),bcInd)=0;
        %put unity on the main diagonal in the positions corresponding to Dirihlet BC
        bc_lin = sub2ind(size(K),bcInd,bcInd);   
        K_bc(bc_lin)=-1;                                        


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculations
        %K_bc = K(bc,bc); 
        %ff_bc = ff(bc);

        %deltaf=ff_bc/norm(Fext(bc)); 
        %u_bc = Regularization(K_bc,ff_bc,SolutionRegType); 
        deltaf=ff_bc/norm(ff); 
        u_bc = Regularization(K_bc,ff_bc,SolutionRegType); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Separation
        % Body1.u(Body1.bc) = Body1.u(Body1.bc) + u_bc(1:Body1.ndof);
        % Body2.u(Body2.bc) = Body2.u(Body2.bc) + u_bc(Body1.ndof + 1:Body1.ndof + Body2.ndof);
        % Body3.u(Body3.bc) = Body3.u(Body3.bc) + u_bc(Body1.ndof + Body2.ndof + 1:end);
        % 
        % Body1.q(Body1.bc) = Body1.q(Body1.bc) + u_bc(1:Body1.ndof);
        % Body2.q(Body2.bc) = Body2.q(Body2.bc) + u_bc(Body1.ndof + 1:Body1.ndof + Body2.ndof);
        % Body3.q(Body3.bc) = Body3.q(Body3.bc) + u_bc(Body1.ndof + Body2.ndof + 1:end);
        Body1.u = Body1.u + u_bc(1:Body1.TotalDofs);
        Body2.u = Body2.u + u_bc(Body1.TotalDofs + 1:Body1.TotalDofs + Body2.TotalDofs);
        Body3.u = Body3.u + u_bc(Body1.TotalDofs + Body2.TotalDofs + 1:end);
        
        Body1.q = Body1.q + u_bc(1:Body1.TotalDofs);
        Body2.q = Body2.q + u_bc(Body1.TotalDofs + 1:Body1.TotalDofs + Body2.TotalDofs);
        Body3.q = Body3.q + u_bc(Body1.TotalDofs + Body2.TotalDofs + 1:end);

        titer=toc;
        titertot=titertot+titer;


        if printStatus(deltaf, u_bc, Re, 1, ii, imax, steps, titertot, Gap)
            break;  
        end 

    end
   


    %Pick nodal displacements from result vector
    xlocName1 = 'xloc' + Body1.ElementType;
    uf1 = Body1.u(Body1.fextInd);
    rf1 = ff(Body1.fextInd);
    rf2 = ff(Body2.fextInd+Body1.TotalDofs);
    rf3 = ff(Body3.fextInd+Body1.TotalDofs+Body2.TotalDofs);

    Results1 = [Results1; Body1.ElementNumber Body1.TotalDofs uf1' rf1'];

    xlocName2 = 'xloc' + Body2.ElementType;
    uf2 = Body2.u(Body2.fextInd); 
    Results2 = [Results2; Body2.ElementNumber Body2.TotalDofs uf2' rf2'];

    xlocName3 = 'xloc' + Body3.ElementType;
    uf3 = Body3.u(Body3.fextInd); 
    Results3 = [Results3; Body3.ElementNumber Body3.TotalDofs uf3' rf3'];
%end    

initial1 = Body1.ForceVectorInit;
initial2 = Body2.ForceVectorInit;
initial3 = Body3.ForceVectorInit;

% % POST PROCESSING ###############################################
hold on
%axis equal
xlabel('\it{X}','FontName','Times New Roman','FontSize',[20])
ylabel('\it{Y}','FontName','Times New Roman','FontSize',[20]),
zlabel('Z [m]','FontName','Times New Roman','FontSize',[20]);
% visualization(Body1,Body1.q,'cyan',true);
% visualization(Body2,Body2.q,'none',true);
% visualization(Body3,Body3.q,'blue',true);
visualization(Body1,Body1.q,'cyan',true);
visualization(Body2,Body2.q,'red',true);
visualization(Body3,Body3.q,'blue',true);

PostProcessing(Body1,Results1,false,false) 
PostProcessing(Body2,Results2,false,false) 
PostProcessing(Body3,Results3,false,false)

CleanTemp(Body1, true)
CleanTemp(Body2, true)
CleanTemp(Body3, true)

