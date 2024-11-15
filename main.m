clc,clear,close all;
format long
ConnectPackages;
% ############ Visualization ##############################################
Ffigplot=0; %  Displacement
% ########## Calculation way ##############################################
sol_acegen = true; % false - Matlab Finite difference, true - AceGen
disp_based = false; % Tensors-based field: displacement (true), position (false)
small_deformation = false; % appoximation theory: infinite small (true), finite (false)
% ########### Element type & data #########################################
Element=3363; % Elements. There are 3243, 3333, 3353, 3363, 34X3 (34103)
ElementData;  
nmesh=1; % Element numbers 
% ########## Problem's data ###############################################
Material=0; % Material models: GOH (0), Neo-Hookean (1), 2- and 5- contant Mooney-Rivlin (2, 5),  K.-S. (3).
Case=1;     % 0 - no load, 1 - elongation load, 2 - bending load 
Area=1;     % 0 - tendon, 1 - rectangular, 2 - circular, 3 - "C" cross-section, 4 - flower 
% ########## Integration way ##############################################
steps = 10; % sub-loading steps
App = 0;    % 0 - standart, 1 and higher - via Green's formula (n is number of appr.)
n_xi = 3;   % number integration point in xi direction
h=10^(-9);  % finite difference scheme step
Pointpic=0; % Picture of integration points in the cross-section
% ########## Problem's data ###############################################
ProblemData; % all information is collected here 
ConnectInnerEnergyFunctions;
% ########### Start the program ###########################################
for n = nmesh % Loop over all defined meshes: nmesh=[1 ....])        
    clear P0f P0 P00 xloc u0 uu ee bc K ff ffc Kc;    % clear previous mesh definitions
    %########## Creates finite element mesh for a simple beam-type structure #############  
    CreateFEMesh;    
    %########## Creates boundary conditions bc #############  
    CreateBC;        
    ndof = sum(bc);      % Number of unconstrained DO  Fs                
    titertot=0;  
    Dvec=[const,H,W,Ln]; % gathering physical and dimensional element data in one vector
    %START NEWTON'S METHOD   
    for i=1:steps
        F_applied=Fmesh(:,i);         
        Re=10^(-4);                   % Stopping criterion for residual
        imax=20;                      % Maximum number of iterations for Newton's method 
        for ii=1:imax    
            tic; 
            [K,ff,Fext] = Kt(q,q0,q0f,u,phim,Phim,F_applied,xloc,nx,nn,nl,MaterialName,Fibers,Dvec,detF0,Gint,Nint,sol_acegen,Element,ElemDofs,DofsAtNode,PosDofs,h);
            Kc = K(bc,bc);            % Eliminate linear constraints from stiffness matrix
            ffc=ff(bc);               % Eliminate linear constraints from force vector
            deltaf=ffc/norm(Fext(bc));% Compute residual
            uc = -Kc\ffc;             % Compute displacements
            u(bc) = u(bc)+uc;         % Add displacement to initial condition
            q(bc) = q(bc)+uc;            
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
        [ux, uy, uz]=Displ(DofsAtNode,u,nn); 
        Case = [Case; n nx ux uy uz];
    end
end
% % POST PROCESSING ###############################################
PostProcessing 