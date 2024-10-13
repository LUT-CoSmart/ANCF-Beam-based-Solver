clc, clear, close all;
disp_based = false;         % what field the tensors are based on: displacemetn (true), position (false)
small_deformation = false;  % appoximation theory: infinite small (true), finite (false)
% ################## Element type & related numbers #######################
Element=3333;                                % Available: 3243, 3333, 3353, 3363, 34X3 (34103)
ElementName = num2str(Element);               % using 'abcd' classification, see in https://doi.org/10.1007/s11071-022-07518-z
Nodes = str2double(ElementName(2));           % Number of nodes            
Dim = str2double(ElementName(end));           % Problem dimensionality     
VecAtNode = str2double(ElementName(3:end-1)); % Vector functions per node  
% ################## Symbolic variables ###################################
syms x y z real;                              % Physical coordinates       
syms xi eta zeta real;                        % Binormalized coordinates   
xi_vec=[xi,eta,zeta].';                       % Binormalized coordinates as a vector 
syms L H W real real;                         % Geometrical variables, element's dimensions 
q = sym('q', [Dim*Nodes*VecAtNode 1], 'real');% Element's position vector of all DoFs
u = sym('u', [Dim*Nodes*VecAtNode 1], 'real');% Element's displacement vector of all DoFs
q0_pos = sym('q0', [Dim*Nodes 1], 'real');    % Element's node position vector in the initial configuration
phi = sym('phi', [Nodes 1], 'real');          % Twist vector of the element's slopes around x axis 
% ################## Basic functions for element ########################## 
% TODO: find a smarter way of derivation "basic functions-element" relation 
BasicFunctions;
switch Element % Find the corresponding basis set to the element
    case 3243,  basis = basis_3243;  required_derivatives = {'x', 'y', 'z'};
    case 3333,  basis = basis_3333;  required_derivatives = {'y', 'z'};      
    case 3353,  basis = basis_3353;  required_derivatives = {'y', 'z', 'yy', 'zz'};
    case 3363,  basis = basis_3363;  required_derivatives = {'y', 'z', 'yz', 'yy', 'zz'};
    case 34103, basis = basis_34103; required_derivatives = {'y', 'z', 'yz', 'yy', 'zz', 'yyz', 'yzz', 'yyy', 'zzz'};
    otherwise   
        error('Unsupported element type!');        
end
ShapeFunctions;
% ################## Element's initial cofiguration #######################
q0 = [];  % Element's DoF vector (all Dofs) in the initial configuration
drdy0 = [0,1,0];
drdz0 = [0,0,1];
for i=1:Nodes 
    A = [1 0 0;
         0 cosd(phi(i)) -sind(phi(i));
         0 sind(phi(i))  cosd(phi(i))];
    drdy = A*drdy0';
    drdz = A*drdz0';
    if ismember('x', required_derivatives)
        drdx=[1;0;0];
    else
        drdx=[];
    end    
    % Collect element's DoF vector (all Dofs) in the initial configuration
    q0 = [q0; q0_pos(Dim*(i-1)+1:Dim*(i-1)+3); drdx; drdy; drdz; ....
          zeros(Dim*VecAtNode-3-length(drdx)-length(drdy)-length(drdz),1)]; % Higher-order terms, their total number = all DoFs - (pos + slopes' lengths)  
end
% ################## Position vectors & tensors ###########################
r0 = Nm_xi*q0;                                       % Position vector in the initial configuration
F0 = jacobian(r0,xi_vec);
if disp_based == true
    uh = Nm_xi*u;                                    % Displacement vector in the actual configuration
    nabla_u = jacobian(uh,xi_vec)*F0^(-1);           % gradient of displacement
    F = eye(3) + nabla_u;                            % deformation gradient via the displacement field
    if small_deformation == true
        dir_name = 'ElementTensors/Displacement/Small';
        E = 0.5*(nabla_u+nabla_u');                  % Green strain tensor based on infinite displacement field  
    else       
        dir_name = 'ElementTensors/Displacement/Finite';
        E = 0.5*(nabla_u+nabla_u'+nabla_u'*nabla_u); % Green strain tensor based on finite displacement field  
    end
    variable = u;
else
    r = Nm_xi*q;                                     % Position vector in the actual configuration
    dir_name = 'ElementTensors/Position'; 
    F = jacobian(r,xi_vec)*F0^(-1);                  % deformation gradient via the position field
    E = 0.5*( F.'*F-eye(3) );                          % Green strain tensor based on position field    
    variable = q;
end
for ii=1:Dim
    for jj=1:Dim
        for kk=1:Dim*Nodes*VecAtNode   
            dEde(ii,jj,kk)=diff(E(ii,jj),variable(kk));
        end
    end
end
% ################## Saving functions #####################################
matlabFunction(Nm_xi, 'file', fullfile('ShapeFunctions', 'Shape_' + string(Element)), 'vars', {L,H,W,xi,eta,zeta});
matlabFunction(F, 'file', fullfile(dir_name,'F_' +  string(Element)), 'vars', {q,u,q0_pos,phi,L,H,W,xi,eta,zeta});
matlabFunction(dEde, 'file', fullfile(dir_name,'dEde_' + string(Element)), 'vars', {q,u,q0_pos,phi,L,H,W,xi,eta,zeta});