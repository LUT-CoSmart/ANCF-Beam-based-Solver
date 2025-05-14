L_Nodes = linspace(0,L,Nodes);
AAA = sym(zeros(length(basis))); 
for i = 0:length(required_derivatives)
    der_Basis = basis;
    if i ~=0
        der_Steps = strlength(required_derivatives{i});          
        for char_idx  = 1:der_Steps 
            symbolic_var = eval(required_derivatives{i}(char_idx));  
            der_Basis =  diff(der_Basis,symbolic_var);
        end       
    end
    for j = 1:Nodes
        k=(j-1)*(length(required_derivatives)+1) + 1;
        AAA(k+i,:) = subs(der_Basis,[x,y,z],[L_Nodes(j), 0, 0]);     
    end
end 
N_x=basis*AAA^-1;      % Shape functions in x coordinatesd
% Mapping x -> xi 
for i=1:length(basis)
    N_xi(i)=simplify(subs(N_x(i),[x,y,z],[(L/2)*(xi+1),(H/2)*eta, (W/2)*zeta]));
    for j = 1:Dim  % Shape function matrix in xi coordinate system
        k=(j-1)+(i-1)*Dim+1;
        Nm_xi(j,k)=N_xi(i);   
        Nm_xi_xi(j,k) = diff(N_xi(i),xi);
        Nm_xi_eta(j,k) = diff(N_xi(i),eta);
        Nm_xi_zeta(j,k) = diff(N_xi(i),zeta);
    end    
end
if Call_shapeFunctions == true
   N_xi   % The presentation of shape functions in isoparametric coordinates 
end    