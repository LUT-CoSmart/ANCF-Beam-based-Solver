% Creates mesh for an intially straight beam structure
addpath('MeshFunctions')
[P0,nloc,phim] = linemesh_ANCF(Element,ElemNodes,n,L,twist_angle,ro);
xloc=xlocAllANCF(ElemNodes,DofsAtNode,nloc);
% get number of nodes (nn), number of elements (nl) and total number of DOFs (nx)
[nn,~] = size(P0);
[nl,~] = size(nloc);
nx = DofsAtNode*nn;          % dofs (no constraints) 
% draw system
u0 = zeros(nx,1);
u = u0;
% create global vector of nodal coordinates
for jj=1:nn
    q((jj-1)*DofsAtNode+1:(jj-1)*DofsAtNode+DofsAtNode)=P0(jj,:); 
end  
% Define initial position
q=q(:);
q0=q;