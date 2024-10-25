function [K_loc,Fe] = ANCFMatlab(MaterialName,K_loc,Fe0,uk,qk,qk0,phik,Phik,Fibers,Dvec,Element,ElemDofs,PosDofs,Gint,Nint,detF0,h)
I_vec=Fe0;              % initialization of the disturbance vector
Fe=Fe_fun(MaterialName,Fe0,uk,qk,qk0,phik,Phik,Fibers,Dvec,Element,ElemDofs,PosDofs,Gint,Nint,detF0); 
for jj = 1:ElemDofs    % cycle over all 
    I_vec(jj)=h;       % disturbance of one coordinate    
    ukh=uk-I_vec;      % displacement vector disturbance
    qkh=qk-I_vec;      % positoin vector disturbance 
    % in the force functions, only one out u and q used, therefore, there is no problem to variate both of them  
    Feh=Fe_fun(MaterialName,Fe0,ukh,qkh,qk0,phik,Phik,Fibers,Dvec,Element,ElemDofs,PosDofs,Gint,Nint,detF0);
    K_loc(:,jj)=(Fe-Feh)/h; 
    I_vec(jj)=0;
end    