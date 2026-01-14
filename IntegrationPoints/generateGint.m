function [Gint,N_int] = generateGint(num, pcirc, wcirc, IntegrationPoints)    

if IntegrationPoints == "Lobatto"
    % [xi,wxi] = Lobatto(num+1);
    [xi,wxi] =  lglnodes(num);
else    
    [xi,wxi] = gauleg2(-1,1,num);
end
% Initialize Gint matrix
Gint = [];
% Loop through all combinations of xiv, etav, and zetav
for i = 1:length(xi)
    for j=1:length(pcirc)
        Gint= [Gint; xi(i), pcirc(j,1), pcirc(j,2), wxi(i) * wcirc(j)];
    end
end
N_int = size(Gint,1); % total number of integration points 
