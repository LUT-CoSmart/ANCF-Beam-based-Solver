function Body = BeamANCFMesh(Body)                                      
    
    % Exctracting info
    DofsAtNode = Body.DofsAtNode; 
    DIM = Body.DIM; 
    ElementNodes = Body.ElementNodes;
    ElementNumber = Body.ElementNumber;

     
    nodes = ElementNumber * (ElementNodes-1) + 1; % number of nodes
    Body.NodeNumber = nodes;
    Body.TotalDofs = DofsAtNode*nodes;          % dofs (no constraints)

    % check fibers
    if Body.Fibers
       Phik=linspace(0,Body.FiberTwist,nodes);
    else
       Phik=zeros(1,nodes); 
    end    

    % check twist
    if isfield(Body, 'Twist')
        phik=linspace(0,Body.Twist.angle,nodes);
        ro = Body.Twist.ro;    
        phi_init = Body.Twist.initial_rot;    
    else
        phik=zeros(1,nodes);
        ro = 0;
        phi_init = 0;
    end    
   
    % Define nodes' positioning 
    xk = linspace(0,Body.Length.X,nodes);    
    yk = ro*cosd(phik+phi_init); % accounting the outer twist
    zk = ro*sind(phik+phi_init); % accounting the outer twist

    % applicable only for beams, used to apply BC locally 
    P00 = [xk' zeros(nodes,1) zeros(nodes,1)]; % simple straight beam-line meshing 
    
    % Slopes' DOFs 
    nullmat = zeros(nodes,3);
    drdy=nullmat;
    drdz=nullmat;
    drdyf=nullmat;
    drdzf=nullmat;

    % rotating slope vectors around x axís
    for i=1:nodes
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

    % Change the beam positioning with dependency on shift and rotation checking shift
    if isfield(Body, 'Shift')
        Shift = [Body.Shift.X, Body.Shift.Y, Body.Shift.Z]; 
    else
        Shift = zeros(1,3);
    end    

    if isfield(Body, 'Rotation')
       alpha = Body.Rotation.X;
       beta = Body.Rotation.Y;
       gamma = Body.Rotation.Z;

       Rx = [1 0 0;
            0 cosd(alpha) -sind(alpha);
            0 sind(alpha) cosd(alpha)];

       Ry = [cosd(beta) 0 sind(beta);
             0 1 0;
            -sind(beta) 0 cosd(beta)];

       Rz = [cosd(gamma) -sind(gamma) 0;
             sind(gamma) cosd(gamma) 0;
             0 0 1];

       R = Rx*Ry*Rz;
    else
       R = eye(3);
    end

    %% TODO: In case of Cosserat add to drdx another vector

    % Here, dependecy is phik = phi/L * xk;
    [Body.Length.Ln,drdx] = SplineLineAlongX(xk,yk,zk,ElementNumber); % Calculating elements' lengths and x-slopes
    if ~Body.Slope_x         
        drdx = []; % It is defined above, but some elements may not use it
    else
        drdx = (R*drdx')';
    end 
    
    drdy = (R*drdy')';
    drdz = (R*drdz')';
    drdyf = (R*drdyf')';
    drdzf = (R*drdzf')';
    
    Positions = (R*[xk' yk' zk']')' + Shift; % update positions due to shifting
    
    % Adding higher-order terms
    HigherOrderDOFs = DofsAtNode-(DIM+size(drdx,2)+DIM+DIM); 
    HigherOrderTerms = zeros(nodes,HigherOrderDOFs);

    % Matrices of all DOFs    
    P0  = [Positions drdx drdy drdz HigherOrderTerms];
    P0f = [Positions drdx drdyf drdzf HigherOrderTerms]; % attention to y- and z-slopes

    % generate element and angles connectivity
    nloc = [];
    phim = [];
    Phim = [];
    for i = 1:ElementNumber
        loc_n = []; % local element's node connectivity
        loc_i = []; % local element's inner angles connectivity
        loc_o = []; % local element's outer angles connectivity
        for j = 1:ElementNodes
            loc_n = [loc_n (i-1)*(ElementNodes-1)+j];
            loc_i = [loc_i phik((i-1)*(ElementNodes-1)+j)];
            loc_o = [loc_o Phik((i-1)*(ElementNodes-1)+j)];
        end
        nloc = [nloc; loc_n];
        phim = [phim; loc_i];
        Phim = [Phim; loc_o];
    end

    Body.nloc = nloc;
    Body.phim = phim;
    Body.Phim = Phim;

    % Creates mesh for an intially straight beam structure
    xlocAllName = "xlocAll" + Body.ElementType;

    Body.xloc=feval(xlocAllName,Body);

    % create global vector of nodal coordinates (q0) and vector for identification of fiber direction in AceGen (q0f) 
    for jj=1:nodes
        q0((jj-1)*DofsAtNode+1:(jj-1)*DofsAtNode+DofsAtNode)=P0(jj,:);
        q0f((jj-1)*DofsAtNode+1:(jj-1)*DofsAtNode+DofsAtNode)=P0f(jj,:);
    end  
    
    % Define initial position
    Body.u = zeros(Body.TotalDofs,1);
    Body.q0=q0(:);
    Body.q=q0(:);
    Body.q0f=q0f(:);

    Rinv = R^(-1);
    nBlock= Body.TotalDofs / size(R,1);   
    Body.RotInv =  kron(eye(nBlock),Rinv);
    Body.Rot = kron(eye(nBlock),R);

    Body.P00 = P00;  % used for application of loadings and BC in local configuration
    