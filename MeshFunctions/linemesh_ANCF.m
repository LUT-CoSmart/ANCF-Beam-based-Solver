function [P0f,P0,nloc,phim,Phim,Ln] = linemesh_ANCF(HigherOrderTerms,Slope_x,ElemNodes,n,L,phi,ro,Phi)                                      
nodes=n*(ElemNodes-1)+1; % number of nodes
if phi ~= 0
    phik=geospace(0,phi,nodes,1)'; % outer twist
else
    phik=zeros(nodes,1); 
end
if Phi ~= 0
    Phik=geospace(0,Phi,nodes,1)';
else
    Phik=zeros(nodes,1);
end    
xk=geospace(0,L,nodes,1)';
% accounting the outer twist
yk=ro*cosd(phik);
zk=ro*sind(phik);
% Basicaly we have dependecy phik = phi/L * xk;
[Ln,drdx] = SplineLineAlongX(xk,yk,zk,n); % Calculating elements' lengths and x-slopes
nullmat = zeros(nodes,3);
drdy=nullmat;
drdz=nullmat;
drdyf=nullmat;
drdzf=nullmat;
% rotating slope vectors
for i=1:nodes
    % form a rotation matrix around x axís
    % counterclockwise rotation 
    A = [1 0 0;
         0 cosd(phik(i)) -sind(phik(i));
         0 sind(phik(i))  cosd(phik(i))];
    drdyk=A*[0;1;0];
    drdzk=A*[0;0;1];
    drdy(i,:)=drdyk';
    drdz(i,:)=drdzk';   
    % Additional rotation for fibers around x axis
    Af = [1 0 0;    
         0 cosd(Phik(i)) -sind(Phik(i));
         0 sind(Phik(i))  cosd(Phik(i))];        
    drdyf(i,:)=(Af*drdyk)';
    drdzf(i,:)=(Af*drdzk)';       
end
HigherOrderTerms_all = zeros(nodes,HigherOrderTerms); % Assumption that dimension is 3                                                      
% Figuring out the usage of x-slope in the chosen element
if Slope_x == false   
   drdx = []; % It is defined above, but some elements may not use it 
end   
P0  = [xk yk zk drdx drdy drdz HigherOrderTerms_all];
P0f = [xk yk zk drdx drdyf drdzf HigherOrderTerms_all]; % attention to y- and z-slopes
% generate element and angles connectivity
nloc = [];
phim = [];
Phim = [];
for i = 1:n
    loc_n = []; % local element's node connectivity
    loc_i = []; % local element's inner angles connectivity
    loc_o = []; % local element's outer angles connectivity
    for j = 1:ElemNodes
        loc_n = [loc_n (i-1)*(ElemNodes-1)+j];
        loc_i = [loc_i phik((i-1)*(ElemNodes-1)+j)];
        loc_o = [loc_o Phik((i-1)*(ElemNodes-1)+j)];
    end
    nloc = [nloc; loc_n];
    phim = [phim; loc_i];
    Phim = [Phim; loc_o];
end