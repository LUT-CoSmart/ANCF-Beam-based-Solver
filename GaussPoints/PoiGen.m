function [pcirc,wcirc]=PoiGen(data,nu2,n)

% Area - points to of the area to deal with;
% n  - how accurate you want to have, the more the better!;
% pic - do you want to see, how points are distributed on the canvas; 
%% ksi-line (doesn't matter wich one)
ksi=@(eta) eta-eta; % all points inside of the cross-section
p = min([n+1,max(cellfun(@numel, data)/2)-1]); % spline degree approximation, chooses between number min numbers of knots-1 (cellfun calculates altogether x and y, therefore, /2)
% and approximation degree +1  (+1 because spapi required at least to be 2)  
%% starting the program
for i=1:nu2
    x=data{i}(:,1);
    y=data{i}(:,2);
    P= [x y];
    m=length(y);
    %% wildely aqdopte parametrization ("cumulative chordal")
    t=zeros(m,1);
    for j=1:m-1
        t(j+1)=t(j)+norm(P(j+1,:)-P(j,:));
    end
    %% approximation 
    S1 = spapi(optknt(t',p),t,x);
    S2 = spapi(optknt(t',p),t,y);
    % premitive of y function 
    % n_i calculation
    if mod(p,2)==1
        n_i=n*p+(p+1)/2;
    else
        n_i=n*p+p/2;
    end
    [tau,w]=lgwt(n_i,-1,1);   % n_i number of Gauss's points
    [tau_n,w_n]=lgwt(n,-1,1); % n number of Gauss's points
    for j=1:m-1
        dt=t(j+1)-t(j);
        for k=1:n_i
            q(j,k)=dt/2*tau(k)+(t(j+1)+t(j))/2;
        end
    end
    point1=[]; % collection of all points
    for j=1:m-1
        dt=t(j+1)-t(j);
        for k=1:n_i
            y_lambda=fnval(S2,q(j,k)); 
            for h=1:n
                hel=fnval(S1,q(j,k)); %help for calculation in point 
                x_lambda=(hel-ksi(y_lambda))/2*tau_n(h)+(hel+ksi(y_lambda))/2;
                S2d=fnval(fnder(S2,1),q(j,k));
                w_lambda=dt/4*w(k)*w_n(h)*(hel-ksi(y_lambda))*S2d;
                point1=[point1; w_lambda x_lambda y_lambda]; 
            end    
        end    
    end
    point{i}=point1;
end
point_all=[];
for i=1:nu2
    point_all=[point_all; point{i}];
end

pcirc(:,2)=point_all(:,2);
pcirc(:,1)=point_all(:,3);
wcirc=point_all(:,1);