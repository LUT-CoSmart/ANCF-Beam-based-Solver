function [K_pp,Fe] = ANCFMatlab(Material,K_pp,Fe0,eek,Dvec,Element,ElemDofs,Gint,Nint,detF0,h)
I_vec=Fe0;              % initialization of the disturbance vector
Fe=Fe_fun(Material,Fe0,eek,Dvec,Element,ElemDofs,Gint,Nint,detF0);   
for jj = 1:ElemDofs    % cycle over all 
    I_vec(jj)=h;       % disturbance of one coordinate    
    eekh=eek-I_vec;    % displacement vector disturbance
    Feh=Fe_fun(Material,Fe0,eekh,Dvec,Element,ElemDofs,Gint,Nint,detF0);   
    K_pp(:,jj)=(Fe-Feh)/h; 
    I_vec(jj)=0;
end    
