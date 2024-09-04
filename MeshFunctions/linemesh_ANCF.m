function [P,nloc] = linemesh_ANCF(Element,ElemNodes,n,L)

% generate nodal coordinates in xy-plane 
nodes=n*(ElemNodes-1)+1; % number of nodes
xk=geospace(0,L,nodes,1)';
yk=zeros(1,nodes)';
zk=zeros(1,nodes)';

nullmat = zeros(nodes,3);
onesvec = ones(nodes,1); 
drdx=nullmat;
drdx(:,1)=onesvec;
drdy=nullmat;
drdy(:,2)=onesvec;
drdz=nullmat;
drdz(:,3)=onesvec;

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
% generate element connectivity
nloc = [];
for i = 1:n
    loc = [];
    for j = 1:ElemNodes
        loc = [loc (i-1)*(ElemNodes-1)+j];
    end
    nloc = [nloc; loc];
end
