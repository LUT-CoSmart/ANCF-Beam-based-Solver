PosDofs = []; % identefication of positional DoFs within the chosen element 
for i = 1:ElemNodes
    PosDofs = [PosDofs (i-1)*DofsAtNode+1:(i-1)*DofsAtNode+3];
end   
if Slope_x 
    length_rx = DIM;
else 
    length_rx = 0;
end    
HigherOrderTerms = DIM*VecAtNode-(DIM+length_rx+DIM+DIM); 
% Creates mesh for an intially straight beam structure
[P0,P0f,nloc,phim,Phim,Ln] = linemesh_ANCF(HigherOrderTerms,Slope_x,ElemNodes,n,L,twist_angle,ro,fiber_twist);
xloc=xlocAllANCF(ElemNodes,DofsAtNode,nloc);
% get number of nodes (nn), number of elements (nl) and total number of DOFs (nx)
[nn,~] = size(P0);
[nl,~] = size(nloc);
nx = DofsAtNode*nn;          % dofs (no constraints) 
% draw system
u = zeros(nx,1);
% create global vector of nodal coordinates
for jj=1:nn
    q((jj-1)*DofsAtNode+1:(jj-1)*DofsAtNode+DofsAtNode)=P0(jj,:); 
    q0f((jj-1)*DofsAtNode+1:(jj-1)*DofsAtNode+DofsAtNode)=P0f(jj,:);
end  
% Define initial position
q=q(:);
q0f=q0f(:);
q0=q;
