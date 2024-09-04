function [Gint,N_int] = generateGint(n,pcirc,wcirc)    
[xiv,wxi]=gauleg2(-1,1,n);
% Initialize Gint matrix
Gint = [];
% Loop through all combinations of xiv, etav, and zetav
for i = 1:n
    for j=1:length(pcirc)
        Gint= [Gint; xiv(i), pcirc(j,1), pcirc(j,2), wxi(i) * wcirc(j)];
    end
end
N_int = size(Gint,1); % total number of integration points 
