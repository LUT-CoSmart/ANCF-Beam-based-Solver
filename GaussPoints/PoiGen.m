function [pcirc,wcirc]=PoiGen(data,nu2,n,pic)
% Area - points to of the area to deal with;
% n  - how accurate you want to have, the more the better!;
% pic - do you want to see, how points are distributed on the canvas; 
%% ksi-line (doesn't matter wich one)
ksi=@(eta) eta-eta; % all points inside of the cross-section
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
    p=2;
    %p=m-1; % degree of "spline"
%     if m-1>2
%         p=3;
%     else
%         p=2;
%     end    
    %p=2; % better for bending and elong ten
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
%% Drawing
if pic==1
    for j=1:length(point_all) 
        xx=point_all(j,2);
        yy=point_all(j,3);
        plot(xx,yy,'.k')
        hold on
        text(xx,yy,num2str(j),'Color','k')
    end
    %% For draw use x1 (not to spoit storage data)
    x1=-1.0:0.1:1.0;
    grid on;
    xlim([-1.1 1.1]);
    ylim([-1.1 1.1]);
    plot(ksi(x1),x1,'-k','LineWidth',2);   
    set(gca,'Fontsize',14);
    set(gca,'FontName','Times New Roman');
end

pcirc(:,2)=point_all(:,2);
pcirc(:,1)=point_all(:,3);
wcirc=point_all(:,1);