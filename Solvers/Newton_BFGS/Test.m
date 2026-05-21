clc,clear,close all;
format long

addpath("MainFunctions");
addpath("MeshFunctions")
addpath("Postprocessing");
addpath('InnerForceFunctions');
addpath(genpath("Solvers"))

Body.Name = "Body";

% ########### Problem data ################################################
Body = DefineElement(Body,"Beam","ANCF",3333,"None"); % 1 - BodyName, 2 - type (beam, plate, etc.), 3 - element name, 4 - modification name (None, EDG, etc.)  
                                                      % ANCF Beam: 3243, 3333, 3343, 3353, 3363, 34X3 (34103) 
Body = Geometry(Body,"Rectangular","Standard","Gauss");  % Cross Sections: Rectangular, Oval, C, Tendon, etc.
                                                      % Integration points of generating line: Gauss, Lobatto  
Body = Materials(Body,"GOH", "optimized_SEE");        % Material models: GOH (GOH, Amir), Neo-Hookean (Neo), 2- and 5- constant Mooney-Rivlin (Mooney2, Mooney5),  Kirhhoff-Saint-Venant (KS).
                                                      % Integration Scheme: Poigen, Standard
                                                     
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
Body.FiniteDiference= "AceGen"; % Calculation of FD: Matlab, Matlab_automatic, AceGen
Body.SolutionBase = "Position"; % Solution-based calculation: Position, Displacement
Body.DeformationType = "Finite"; % Deformation type: Finite, Small
Body = AddTensors(Body);

%% TODO: rebuild CreateMex, it addresses the wrong folder
Body.mex = false;

% ########## Boundary Conditions ##########################################
% Force 
Force.Maginutude.X = 2e6;  % Elongation

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
 
% % %####################### Solving ######################################## 
steps = 300;  % sub-loading steps
titertot=0;  
Re=10^(-4);           % Stopping criterion for residual
imax=20;              % Maximum number of iterations for Newton's method 

%START NEWTON'S METHOD   

Body = SubLoading(Body, 1, steps, "linear"); 

myFx = @(x)myfun(Body,x);
x0 = zeros(Body.ndof,1);

maxIter = 100;

m = 20;
re = 1e-5;

n = length(x0);
Sm = zeros(n,m);
Ym = zeros(n,m);
[f0,g0]=feval(myFx,x0)

alpha = ones(n,1);
x1 = x0 - g0* alpha


k = 1;
while (k < maxIter) && norm(g0)>re

    fnorm = norm(g0)   

    s0 = x1-x0;
    y0 = g1-g0;

    hdiag = s0'*y0/(y0'*y0);
    p = zeros(length(g0),1);
    if (k<=m)
        % update S,Y
        Sm(:,k) = s0;
        Ym(:,k) = y0;
        % never forget the minus sign
        p = -getHg_lbfgs(g1,Sm(:,1:k),Ym(:,1:k),hdiag); 
    elseif (k>m)
        Sm(:,1:(m-1))=Sm(:,2:m);
        Ym(:,1:(m-1))=Ym(:,2:m);
        Sm(:,m) = s0;
        Ym(:,m) = y0;    
        % never forget the minus sign
        p = -getHg_lbfgs(g1,Sm,Ym,hdiag);
    end  

    % line search
    [alpha,fs,gs]= strongwolfe(myFx,p,x1,f1,g1);
    x0 = x1;
    g0 = g1;
    x1 = x1 + alpha*p;
    f1 = fs;
    g1 = gs;
    k =k + 1;  
end
% 
% % 
% % Body.u(Body.bc) = Body.u(Body.bc)+u_bc;         % Add displacement to previous one
% % Body.q(Body.bc) = Body.q(Body.bc)+u_bc;         % change the global positions
% % 
% % titer=toc;
% % titertot=titertot+titer;   
% 
