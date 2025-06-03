clc, clear, close all;
Call_shapeFunctions = false;% To create AceGen-generated functions, it might be necessary to demonstrate shape functinos 
write_files = true;         % Do we need to write any file?
disp_based = true;          % What field the tensors are based on: displacement (true), position (false)
small_deformation = false;  % Appoximation theory: infinite small (true), finite (false)
% ################## Element type & related numbers #######################
Element=3333;                                 % Available: 3243, 3333, 3343, 3353, 3363, 34X3 (34103)
ElementName = num2str(Element);               % using 'abcd' classification, see in https://doi.org/10.1007/s11071-022-07518-z
Nodes = str2double(ElementName(2));           % Number of nodes            
Dim = str2double(ElementName(end));           % Problem dimensionality     
VecAtNode = str2double(ElementName(3:end-1)); % Vector functions per node  
% ################## Symbolic variables ###################################
syms x y z real;                              % Physical coordinates       
syms xi eta zeta real;                        % Binormalized coordinates   
xi_vec=[xi;eta;zeta];                         % Binormalized coordinates as a vector 
syms L H W real real;                         % Geometrical variables, element's dimensions 
q = sym('q', [Dim*Nodes*VecAtNode 1], 'real');% Element's position vector of all DoFs
u = sym('u', [Dim*Nodes*VecAtNode 1], 'real');% Element's displacement vector of all DoFs
q0pos = sym('q0pos', [Dim*Nodes 1], 'real');  % Element's node position vector in the initial configuration
q0x = sym('q0x', [Dim*Nodes 1], 'real');      % Element's coordinate along x direction, in case we have an axial curve
phi = sym('phi', [Nodes 1], 'real');          % Twist vector of the element's slopes around x axis 
a0 = sym('a0', [Dim 1],'real');               % Fibres' direction vector
Phi = sym('Phi', [Nodes 1], 'real');          % Pre-twist of fiber vector around x axis of element's slopes 
F_brick_inv = [2/L 0 0; 0 2/H 0; 0 0 2/W];    % Gradient of brick (dF_00/d_xi)^-1 used for fibers adjustments
% ################## Basic functions for element ########################## 
BasicFunctions;
switch Element % Find the corresponding basis set to the element
    case 3243,  basis = basis_3243;  required_derivatives = {'x', 'y', 'z'};
    case 3333,  basis = basis_3333;  required_derivatives = {'y', 'z'};      
    case 3343,  basis = basis_3343;  required_derivatives = {'y', 'z', 'yz'};              
    case 3353,  basis = basis_3353;  required_derivatives = {'y', 'z', 'yy', 'zz'};
    case 3363,  basis = basis_3363;  required_derivatives = {'y', 'z', 'yz', 'yy', 'zz'};
    case 34103, basis = basis_34103; required_derivatives = {'y', 'z', 'yz', 'yy', 'zz', 'yyz', 'yzz', 'yyy', 'zzz'};
    otherwise   
        error('Unsupported element type!');        
end
ShapeFunctions;
% ################## Element's initial cofiguration #######################
q0 = [];   % Element's DoF vector (all Dofs) in the initial configuration
q0f = []; % Element's DoF vector (all Dofs) in the initial fibers' configuration for the length reshaping
q0posAndrx = [];
drdy_brick = [0;1;0];
drdz_brick = [0;0;1];
for i=1:Nodes 
    A = [1 0 0;
         0 cosd(phi(i)) -sind(phi(i));
         0 sind(phi(i))  cosd(phi(i))];
    Af = [1 0 0;    % Additional rotation for fibers around x axis 
         0 cosd(Phi(i)) -sind(Phi(i));
         0 sind(Phi(i))  cosd(Phi(i))];    
    drdy_curve = A*drdy_brick;
    drdz_curve = A*drdz_brick;
    if ismember('x', required_derivatives)
        drdx_curve= [q0x(Dim*(i-1)+1:Dim*(i-1)+3)];
    else
        drdx_curve=[];
    end    
    % Collect element's DoF vector (all Dofs) in the initial configuration
    HigherOrderTerms = zeros(Dim*VecAtNode-(Dim+length(drdx_curve)+Dim+Dim ),1); % Number = Dim * VecAtNodes - (pos + slopes' lengths) 
    q0 = [q0; q0pos(Dim*(i-1)+1:Dim*(i-1)+Dim); drdx_curve; drdy_curve; drdz_curve; HigherOrderTerms];
    q0f = [q0f; q0pos(Dim*(i-1)+1:Dim*(i-1)+Dim); drdx_curve; Af*drdy_curve; Af*drdz_curve; HigherOrderTerms];    

    q0posAndrx = [q0posAndrx; q0pos(Dim*(i-1)+1:Dim*(i-1)+Dim);  drdx_curve];    
end

% % ################## Position vectors & tensors ###########################
I = eye(3);                                          % identity tensor
r0 = Nm_xi*q0;                                       % Position vector in the initial configuration
F0 = jacobian(r0,xi_vec);
if disp_based == true
    uh = Nm_xi*u;                                    % Displacement vector in the actual configuration
    nabla_u = jacobian(uh,xi_vec)*F0^(-1); % gradient of displacement
    F = I + nabla_u;                            % deformation gradient via the displacement field
    if small_deformation == true
        dir_name = 'Displacement/Small';
        E = 0.5*(nabla_u+nabla_u');                  % Green strain tensor based on infinite displacement field  
    else       
        dir_name = 'Displacement/Finite';
        E = 0.5*(nabla_u+nabla_u'+nabla_u'*nabla_u); % Green strain tensor based on finite displacement field  
    end
    variable = u;
else
    r = Nm_xi*q;                                     % Position vector in the actual configuration
    dir_name = 'Position'; 
    F = jacobian(r,xi_vec)*F0^(-1);                  % deformation gradient via the position field
    E = 0.5*( F'*F - I );                            % Green strain tensor based on position field    
    variable = q;
end
for ii=1:Dim
    for jj=1:Dim
        for kk=1:Dim*Nodes*VecAtNode   
            dEde(ii,jj,kk)=diff(E(ii,jj),variable(kk));
        end
    end
end

% ################## Function for fibers' transformation ##################
% Here is assumed that fibers are measured in a straght brick sample, 
% which makes the identification of their direction easier.
% Another assumption: the element is firstly pre-deformed,
% then the fibers are pre-twisted around x axis in the counterclockwise direction. 
r0f = Nm_xi*q0f;                       % Position vector in the initial configuration
F0f = jacobian(r0f,xi_vec);            % Gradient of "fibers" in the pre-deformed configuration  
F_fibers = simplify(F0f * F_brick_inv);% Total fiber reshape from rectangular brick to the intial pre-deformed configuration 
a0_fib=F_fibers*a0;                    % Reshaping the fibers
a0_fib=a0_fib/norm(a0_fib);            % Normalization, such as the length can change 

% ################## Saving functions #####################################
if write_files == true
    matlabFunction(Nm_xi, 'file', fullfile('ShapeFunctions','Shape_' + string(Element)), 'vars', {L,H,W,xi,eta,zeta});
    matlabFunction(Nm_xi_xi, 'file', fullfile('ShapeFunctions','Shape_xi_' + string(Element)), 'vars', {L,H,W,xi,eta,zeta});
    matlabFunction(Nm_xi_eta, 'file', fullfile('ShapeFunctions','Shape_eta_' + string(Element)), 'vars', {L,H,W,xi,eta,zeta});
    matlabFunction(Nm_xi_zeta, 'file', fullfile('ShapeFunctions','Shape_zeta_' + string(Element)), 'vars', {L,H,W,xi,eta,zeta});
    
    matlabFunction(a0_fib, 'file', fullfile('FiberVectors','a0_fib_' + string(Element)), 'vars', {a0,q0posAndrx,phi,Phi,L,H,W,xi,eta,zeta});
    matlabFunction(F, 'file', fullfile(dir_name,'F_' +  string(Element)), 'vars', {q,u,q0posAndrx,phi,L,H,W,xi,eta,zeta});
    matlabFunction(dEde, 'file', fullfile(dir_name,'dEde_' + string(Element)), 'vars', {q,u,q0posAndrx,phi,L,H,W,xi,eta,zeta});
end
