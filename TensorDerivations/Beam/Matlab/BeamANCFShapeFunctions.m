clc, clear, close all;
Call_shapeFunctions = true;   % To create AceGen-generated functions, it might be necessary to demonstrate shape functinos 
write_files = true;            % Do we need to write any file?

% ################## Element type & related numbers #######################
Element=3333;                                 % Available: 3243, 3333, 3343, 3353, 3363, 34103 (34X3)
ElementName = num2str(Element);               % using 'abcd' classification, see in https://doi.org/10.1007/s11071-022-07518-z
Nodes = str2double(ElementName(2));           % Number of nodes            
Dim = str2double(ElementName(end));           % Problem dimensionality     
VecAtNode = str2double(ElementName(3:end-1)); % Vector functions per node  
NDoF = Dim*VecAtNode*Nodes;
% ################## Symbolic variables ###################################
syms x y z real;                                  % Physical coordinates       
syms xi eta zeta real;                            % Binormalized coordinates   
xi_vec=[xi;eta;zeta];                             % Binormalized coordinates as a vector 
syms L H W real;                             % Geometrical variables, element's dimensions 

% ################## Basic functions for element ########################## 
basis_3243 =[1, x, z, x^2, y, x*y, x*z, x^3];
basis_3333 =[1, x, z, x^2, y, x*y, x*z, x^2*y, x^2*z]; 
basis_3343 =[1, x, z, x^2, y, x*y, x*z, x^2*y, x^2*z, y*z, x*y*z, x^2*y*z];
basis_3353 =[1, x, z, x^2, y, x*y, x*z, x^2*y, x^2*z, y^2, x*y^2, x^2*y^2, z^2, x*z^2, x^2*z^2];
basis_3363 =[1, x, z, x^2, y, x*y, x*z, x^2*y, x^2*z, y^2, x*y^2, x^2*y^2, z^2, x*z^2, x^2*z^2, y*z, x*y*z, x^2*y*z];
basis_34103=[1, x, z, x^2, y, x*y, x*z, x^2*y^2, x^2*y^3, x^3*y^2, x^3*y^3, x^2*z^2, x^2*z^3, x^3*z^2, x^3*z^3, y*z, x*y^2, x^2*y, x*y^3, x^3*y, x*z^2, x^2*z, x*z^3, x^3*z, y*z^2, y^2*z, x^3, y^2, y^3, z^2, z^3, x*y*z^2, x*y^2*z, x^2*y*z, x^3*y*z, x^2*y*z^2, x^2*y^2*z, x^3*y*z^2, x^3*y^2*z, x*y*z];

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

% ################## Element's Shape Functions #######################
L_Nodes = linspace(0,L,Nodes);
AAA = sym(zeros(length(basis))); 
for i = 0:length(required_derivatives)
    der_Basis = basis;
    if i ~=0
        der_Steps = strlength(required_derivatives{i});          
        for char_idx  = 1:der_Steps 
            symbolic_var = eval(required_derivatives{i}(char_idx));  
            der_Basis =  diff(der_Basis,symbolic_var);
        end       
    end
    for j = 1:Nodes
        k=(j-1)*(length(required_derivatives)+1) + 1;
        AAA(k+i,:) = subs(der_Basis,[x,y,z],[L_Nodes(j), 0, 0]);     
    end
end 
N_x=basis*AAA^-1;      % Shape functions in x coordinatesd
% Mapping x -> xi 
for i=1:length(basis)
    N_xi(i)=simplify(subs(N_x(i),[x,y,z],[(L/2)*(xi+1),(H/2)*eta, (W/2)*zeta]));
    for j = 1:Dim  % Shape function matrix in xi coordinate system
        k=(j-1)+(i-1)*Dim+1;
        Nm_xi(j,k)=N_xi(i);   
        Nm_xi_xi(j,k) = diff(N_xi(i),xi);
        Nm_xi_eta(j,k) = diff(N_xi(i),eta);
        Nm_xi_zeta(j,k) = diff(N_xi(i),zeta);
    end    
end
if Call_shapeFunctions == true
   N_x
   N_xi   % The presentation of shape functions in isoparametric coordinates 
end    

% ################## Saving functions #####################################
if write_files == true
    
    dir_name = "ShapeFunctions/" + ElementName;

    if ~isfolder(dir_name)
        mkdir(dir_name);
    end

    matlabFunction(Nm_xi, 'file', fullfile(dir_name,'Shape_'), 'vars', {L,H,W,xi,eta,zeta});
    matlabFunction(Nm_xi_xi, 'file', fullfile(dir_name,'Shape_xi_'), 'vars', {L,H,W,xi,eta,zeta});
    matlabFunction(Nm_xi_eta, 'file', fullfile(dir_name,'Shape_eta_'), 'vars', {L,H,W,xi,eta,zeta});
    matlabFunction(Nm_xi_zeta, 'file', fullfile(dir_name,'Shape_zeta_'), 'vars', {L,H,W,xi,eta,zeta});
  
end
