clc, clear, close all;

import casadi.*

ElementName = 'ANCF3333';
Material = 'TOM';

if strcmp(Material, 'KS') || strcmp(Material, 'Neo')
    Dvec = SX.sym('p', 5, 1); % Vector of parameters of KS: [E, nu, H, W, L]
elseif strcmp(Material, 'GOH') 
    Dvec = SX.sym('p', 11, 1); % Vector of parameters of GOH: [c10, k1, k2, kappa, a0, d, H, W, L]    
elseif strcmp(Material, 'TOM') 
    Dvec = SX.sym('p', 12, 1);  
else    

end

% Identity tensor
DIM = 3;
I = SX.eye(DIM);

% Nodal coordinates (actual)
q = SX.sym('q', 27, 1);

% Nodal coordinates (initial)
q0 = SX.sym('q0', 27, 1);

% Nodal coordinates of pre-twisted fibers (initial)
q0f = SX.sym('q0f', 27, 1);

% Local coordinates
xi   = SX.sym('xi',1,1);
eta  = SX.sym('eta',1,1);
zeta = SX.sym('zeta',1,1);
xiv = [xi; eta; zeta];

% Geometrical parameters

H = Dvec(end-2);
W = Dvec(end-1);
L = Dvec(end);

% Shape functions 
N_xi = SX.zeros(1,9);
% 1 Node
N_xi(1) = xi*(xi - 1)/2;
N_xi(2) = H*eta*xi*(xi - 1)/4;
N_xi(3) = W*xi*zeta*(xi - 1)/4;
% 2 Node
N_xi(4) = 1 - xi^2;
N_xi(5) = -H*eta*(xi^2 - 1)/2;
N_xi(6) = -W*zeta*(xi^2 - 1)/2;
% 3 Node
N_xi(7) = xi*(xi + 1)/2;
N_xi(8) = H*eta*xi*(xi + 1)/4;
N_xi(9) = W*xi*zeta*(xi + 1)/4;

% Shape function matrix
Smxi = [N_xi(1)*I N_xi(2)*I N_xi(3)*I ...
        N_xi(4)*I N_xi(5)*I N_xi(6)*I ...
        N_xi(7)*I N_xi(8)*I N_xi(9)*I];

% Position vectors
r0 = Smxi*q0;
r0f = Smxi*q0f; % pretwist of fibers from measured rectangular configuration
r  = Smxi*q;
% Deformation gradient
F = jacobian(r,xiv)*(jacobian(r0,xiv))^(-1);
F0f = jacobian(r0f,xiv) * [2/L 0 0; 0 2/H 0; 0 0 2/W]; % Gradient of brick (dF_00/d_xi)^-1 used for fibers adjustments 
% Green-Lagrange strain
E = 0.5*(F'*F-I);

if strcmp(Material, 'KS')
    Young = Dvec(1);
    Poisson = Dvec(2);
    
    lambda = Young*Poisson/((1+Poisson)*(1-2*Poisson));
    mu =  Young/(2*(1+Poisson));
    Psi = 0.5*lambda*trace(E)^2 + mu*sumsqr(E);

elseif strcmp(Material, 'Neo')
    % Material  parameters
    mu = Dvec(1);
    d  = Dvec(2);    
    I1  = trace(F'*F);
    J  = det(F);
    Psi = mu/2 * (I1 - 3) + 1/d*(J-1)^2;
elseif strcmp(Material, 'GOH')   
    c10 = Dvec(1);
    k1 = Dvec(2);
    k2 = Dvec(3);
    kappa = Dvec(4);
    a0 = Dvec(5:7);     % fibers measured in a straght brick sample, 
    a0_fib = F0f*a0(:); % adjustments in the possible pre-twist 
    d = Dvec(8);
    
    C = F'*F;
    I1  = trace(C);
    J  = det(F);
    I4 = a0_fib' * C * a0_fib;
    I4 = kappa*I1 + (1 - 3*kappa) * I4 - 1; 
    Psi = c10/2 * (I1 - 3) + k1/(2*k2) * (exp(k2*I4^2) - 1) + 1/d*(J-1)^2; 
elseif  strcmp(Material, 'TOM')   
    muphi = Dvec(1);
    Ephi  = Dvec(2);
    a     = Dvec(3);
    b     = Dvec(4);
    c     = Dvec(5);
    
    a0     = Dvec(6:8);
    a0_fib = F0f*a0(:);
    
    d = Dvec(9);
    
    C  = F.' * F;
    I1 = trace(C);
    J  = det(F);
    I4 = a0_fib.' * C * a0_fib;
    
    
    %% Denominators
    denom1 = (b-a)*(c-a);
    denom2 = (c-a)*(b-c);
    denom3 = (b-a)*(b-c);
    
    
    %% A
    A1 = -a^2 / denom1;
    A2 =  c^2 / denom2 + A1;
    A3 = -1;
    
    
    %% B
    B1 = 2*a*log(a) / denom1;
    B2 = B1 - 2*c*log(c) / denom2;
    B3 = B2 + 2*b*log(b) / denom3;
    
    
    %% C
    C1 =  1 / denom1;
    C2 = -1 / denom3;
    C3 =  0;
    
    
    %% D
    D1 = -2*a / denom1;
    D2 =  2*b / denom3;
    D3 =  0;
    
    
    %% G
    G1 = (a^2*log(a) - 2.5*a^2) / denom1;
    
    G2 = G1 ...
       - c^2*log(c) / denom2 ...
       + 2.5*c^2 / denom2;
    
    G3 = G2 ...
       + b^2*log(b) / denom3 ...
       - 2.5*b^2 / denom3;
    
    Phi1 = A1/2 * log(I4) ...
         + (B1-D1) * sqrt(I4) ...
         + C1/2 * I4 ...
         + D1/2 * sqrt(I4)*log(I4) ...
         + G1;
    
    Phi2 = A2/2 * log(I4) ...
         + (B2-D2) * sqrt(I4) ...
         + C2/2 * I4 ...
         + D2/2 * sqrt(I4)*log(I4) ...
         + G2;
    
    Phi3 = A3/2 * log(I4) ...
         + (B3-D3) * sqrt(I4) ...
         + C3/2 * I4 ...
         + D3/2 * sqrt(I4)*log(I4) ...
         + G3;
    
    
    %% CasADi piecewise function
    Psi_fiber = if_else(I4 < a^2, ...
                        0, ...
                 if_else(I4 <= c^2, ...
                         Ephi*Phi1, ...
                 if_else(I4 <= b^2, ...
                         Ephi*Phi2, ...
                         Ephi*Phi3)));
    
    
    %% Total energy
    Psi = muphi/2 * (I1 - 3) ...
        + Psi_fiber ...
        + 1/d * (J - 1)^2;
else    

end

Ke = hessian(Psi,q);
Fe = gradient(Psi,q);

% Saving
FunctionName = [ElementName Material];
ElementFun = Function(FunctionName, {q,q0,q0f,Dvec,xi,eta,zeta}, {Ke,Fe},{'q','q0','q0f','p','xi','eta','zeta'},{'Ke','Fe'});
SaveName = [FunctionName '.casadi'];
ElementFun.save(SaveName);

%% TODO: consider this option after adding mex generator to the code
% options = struct('mex',true);
% ElementFun.generate('Element54_codegen.c',options);
% mex Element54_codegen.c -largeArrayDims