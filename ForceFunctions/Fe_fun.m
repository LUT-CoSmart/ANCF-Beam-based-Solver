function Fe = Fe_fun(Material,Fe0,ee,Dvec,Element,ElemDofs,Gint,Nint,detF0) 
L = Dvec(end);   % element length
Fe = Fe0;
for ii=1:Nint    % integration all over the element's volume
    Fe = Fe + FedV(Material,Fe0,ee,Dvec,Element,ElemDofs,Gint(ii,1),Gint(ii,2),Gint(ii,3))*Gint(ii,4);
end

Fe=Fe*1/2*L*detF0;
