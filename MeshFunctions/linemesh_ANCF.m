function [P,nloc,phim] = linemesh_ANCF(Element,ElemNodes,n,L,phi,ro)
                                       % generate nodal coordinates in xy-plane 
nodes=n*(ElemNodes-1)+1; % number of nodes
if phi ~= 0
    phik=geospace(0,phi,nodes,1)'; % outer twist
else
    phik=zeros(nodes,1);
end

xk=geospace(0,L,nodes,1)';
% accounting the outer twist
yk=ro*cosd(phik);
zk=ro*sind(phik);

nullmat = zeros(nodes,3);
onesvec = ones(nodes,1); 
drdx=nullmat;
drdx(:,1)=onesvec;
drdy=nullmat;
drdz=nullmat;
% rotating slope vectors
for i=1:nodes
    % form a rotation matrix around x axís
    A = [1 0 0;
         0 cosd(phik(i)) -sind(phik(i));
         0 sind(phik(i))  cosd(phik(i))];
    % counterclockwise rotation 
    drdyk=A*[0 1 0]';
    drdzk=A*[0 0 1]';
    drdy(i,:)=drdyk.';
    drdz(i,:)=drdzk.';
end

%% TODO: implemnt like it is done in tensor derivation function (more general w/o specification here by elements) 
if Element==3243
   P = [xk yk zk drdx drdy drdz];
elseif Element==3333
   P = [xk yk zk drdy drdz];
elseif Element==3353
   P = [xk yk zk drdy drdz nullmat nullmat];
elseif Element==3363
   P = [xk yk zk drdy drdz nullmat nullmat nullmat];
elseif Element==34103
   P = [xk yk zk drdy drdz nullmat nullmat nullmat nullmat nullmat nullmat nullmat];
end

% generate element and angles connectivity
nloc = [];
phim = [];
for i = 1:n
    loc_n = []; % local element's node connectivity
    loc_p = []; % local element's angles connectivity
    for j = 1:ElemNodes
        loc_n = [loc_n (i-1)*(ElemNodes-1)+j];
        loc_p = [loc_p phik((i-1)*(ElemNodes-1)+j)];
    end
    nloc = [nloc; loc_n];
    phim = [phim; loc_p];
end
