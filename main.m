clc,clear,close all;
format long
addpath('SupportFunctions');
% ############ Visualization ##############################################
Ffigplot=0; %  Displacement
% ########### Element type & data #########################################
Element=3333; % Elements. There are 3243, 3333, 3353, 3363, 34X3 (34103)
ElementData;  
nmesh=1; % Element numbers 
% ########## Problem's data ###############################################
Material=1; % Materials. There are Neo-Hookean (1), 2 - and 5 - contants Mooney-Rivlin (2, 5), GOH (0), K.-S. (3).
Case=1;     % 0 - no load, 1 - elongation load, 2 - bending load 
Area=1;     % 0 - tendon, 1 - rectangular, 2 - circular, 3 - "C" cross-section, 4 - flower 
BC=0;       % 0 - reduced, 1 - full 
% ########## Integration way ##############################################
steps = 10; % sub-loading steps
App = 0;    % 0 - standart, 1 and higher - via Green's formula (n is number of appr.)
n_xi = 3;   % number integration point in xi direction
sol = 1;    % 0 - standart, 1- AceGen
h=10^(-9);  % finite difference scheme step
Pointpic=0; % Picture of integration points in the cross-section
% ########## Problem's data ###############################################
ProblemData; % all information is collected here  
% ########### Start the program ###########################################
for n = nmesh % Loop over all defined meshes: nmesh=[1 ....])     
    clear P P0 xloc u0 uu ee bc f1 K ff ffcs Kc Kcs;    % clear previous mesh definitions
    %########## Creates finite element mesh for a simple beam-type structure #############  
    CreateFEMesh;
    %########## Creates boundary conditions bc #############  
    CreateBC;        
    ndof = sum(bc);     % Number of unconstrained DO  Fs                
    titertot=0;    
    %START NEWTON'S METHOD   
    for i=1:length(Fmesh(1,:))
        F_applied=Fmesh(:,i);         
        Re=10^(-4);           % Stopping criterion for residual
        imax=20;              % Maximum number of iterations for Newton's method 
        for ii=1:imax    
            tic; 
            [K,ff,Fext] = Kt(ee,ee0,F_applied,xloc,nx,nn,nl,Material,Dvec,detF0,Gint,Nint,sol,Element,ElemDofs,DofsAtNode,h);
            %Eliminate linear constraints
            Kc = K(bc,bc); 
            ffc=ff(bc);         
            deltaf=ffc/norm(Fext(bc));  % compute residual
            uuc = -Kc\ffc;          % compute displacements
            uu(bc) = uu(bc)+uuc;    % add displacement to initial condition
            ee(bc) = ee(bc)+uuc;
            
            titer=toc;
            titertot=titertot+titer;   
            if abs(deltaf)<Re*ones(ndof,1)
                disp(['Solution is found by ' num2str(ii) ' iterations. Total CPU-time: ' num2str(titertot)])
                break
            elseif ii==imax 
                disp(['The solution is not found. The maximum number of iterations is reached. Total CPU-time: ' num2str(titertot)])
            else     
                disp('wait...')
            end              
        end           
        %Pick nodal displacements from result vector
        [ux, uy, uz]=Displ(DofsAtNode,uu,nn); 
        Case = [Case; n nx ux uy uz];
    end
end
% % POST PROCESSING ###############################################
PostProcessing 