clc, clear, close all;
Call_shapeFunctions = false;   % To create AceGen-generated functions, it might be necessary to demonstrate shape functinos 
write_files = true;            % Do we need to write any file?

% ################## Element type & related numbers #######################
Element=3363;                                 % Available: 3243, 3333, 3343, 3353, 3363, 34X3 (34103)
ElementName = num2str(Element);               % using 'abcd' classification, see in https://doi.org/10.1007/s11071-022-07518-z
Nodes = str2double(ElementName(2));           % Number of nodes            
Dim = str2double(ElementName(end));           % Problem dimensionality     
VecAtNode = str2double(ElementName(3:end-1)); % Vector functions per node  
NDoF = Dim*VecAtNode*Nodes;
% ################## Symbolic variables ###################################
syms x y z real;                                  % Physical coordinates       
syms xi eta zeta real;                            % Binormalized coordinates   
xi_vec=[xi;eta;zeta];                             % Binormalized coordinates as a vector 
syms L H W real real;                             % Geometrical variables, element's dimensions 
q = sym('q', [Dim*Nodes*VecAtNode 1], 'real');    % Element's position vector of all DoFs
u = sym('u', [Dim*Nodes*VecAtNode 1], 'real');    % Element's displacement vector of all DoFs
q0 = sym('q0', [Dim*Nodes*VecAtNode 1], 'real');  % Element's position vector of all DoFs (initial configuration)
a0 = sym('a0', [Dim 1],'real');                   % Fibres' direction vector
Phi = sym('Phi', [Nodes 1], 'real');              % Pre-twist of fiber vector around x axis of element's slopes 
F_brick_inv = [2/L 0 0; 0 2/H 0; 0 0 2/W];        % Gradient of brick (dF_00/d_xi)^-1 used for fibers adjustments

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
q0f = q0; % Element's DoF vector (all Dofs) in the initial fibers' configuration for the length reshaping
q0s = q0; % Substituted DoFs with zeros for high-ordered terms, assuming the configuration is not pre-stressed

if ismember('x', required_derivatives)
   drdx_len = 3;
else
   drdx_len = 0;
end  

% Number = Dim * VecAtNodes - (pos + slopes' lengths) 
HigherOrderTerms = zeros(Dim*VecAtNode-(Dim+drdx_len+Dim+Dim),1); 

for i=1:Nodes 

    skip_dofs = Dim*VecAtNode*(i-1) + Dim + drdx_len;
    drdy = q0( skip_dofs + 1:skip_dofs+3); 
    drdz = q0( skip_dofs + 4:skip_dofs+6);

    if ~isempty(HigherOrderTerms)  % substitution with zeros to speed up funciton derivations 
        q0f(skip_dofs+7:Dim*VecAtNode*i) = HigherOrderTerms;
        q0s(skip_dofs+7:Dim*VecAtNode*i) = HigherOrderTerms; 
    end

    Af = [1 0 0;    % Additional rotation for fibers around x axis 
         0 cosd(Phi(i)) -sind(Phi(i));
         0 sind(Phi(i))  cosd(Phi(i))];    

    q0f( skip_dofs + 1:skip_dofs+6) = [Af*drdy;Af*drdz]; 
end

% ################## Position vectors & tensors ###########################
dir_name = "ElementFunctions/" + ElementName;
r0 = Nm_xi*q0s;                            % Position vector in the initial configuration
r = Nm_xi*q;                               % Position vector in the actual configuration
uh = Nm_xi*u;                              % Displacement vector in the actual configuration
F0 = jacobian(r0,xi_vec);

F_xi = jacobian(r,xi_vec);          % deformation gradient in local coordinates
nabla_u_xi = jacobian(uh,xi_vec);          % Gradient of displacement in local coordinates

% ################## Function for fibers' transformation ##################
% Here is assumed that fibers are measured in a straght brick sample, 
% which makes the identification of their direction easier.
% element is pre-deformed, then the fibers are pre-twisted around x axis in the counterclockwise direction. 
r0f = Nm_xi*q0f;                       % Position vector in the initial configuration
F0f = jacobian(r0f,xi_vec);            % Gradient of "fibers" in the pre-deformed configuration  
F_fibers = simplify(F0f * F_brick_inv);% Total fiber reshape from rectangular brick to the intial pre-deformed configuration 
a0_fib=F_fibers*a0;                    % Reshaping the fibers
a0_fib=a0_fib/norm(a0_fib);            % Normalization, such as the length can change 

% ################## Saving functions #####################################
if write_files == true
    
    dir_name_2 = "ShapeFunctions/" + ElementName;

    if ~isfolder(dir_name_2)
        mkdir(dir_name_2);
    end

    matlabFunction(Nm_xi, 'file', fullfile(dir_name_2,'Shape_'), 'vars', {L,H,W,xi,eta,zeta});
    matlabFunction(Nm_xi_xi, 'file', fullfile(dir_name_2,'Shape_xi_'), 'vars', {L,H,W,xi,eta,zeta});
    matlabFunction(Nm_xi_eta, 'file', fullfile(dir_name_2,'Shape_eta_'), 'vars', {L,H,W,xi,eta,zeta});
    matlabFunction(Nm_xi_zeta, 'file', fullfile(dir_name_2,'Shape_zeta_'), 'vars', {L,H,W,xi,eta,zeta});

    if ~isfolder(dir_name)
        mkdir(dir_name);
    end

    matlabFunction(a0_fib, 'file', fullfile(dir_name,'a0_fib'), 'vars', {a0,q0,Phi,L,H,W,xi,eta,zeta});
    matlabFunction(F_xi, 'file', fullfile(dir_name,'F_xi'), 'vars', {q,L,H,W,xi,eta,zeta});    
    matlabFunction(F0, 'file', fullfile(dir_name,'F0'), 'vars', {q0,L,H,W,xi,eta,zeta});
    matlabFunction(nabla_u_xi , 'file', fullfile(dir_name,'nabla_u_xi'), 'vars', {u,L,H,W,xi,eta,zeta});
end
