function incenter=getQuadInCenter(pts)
% determine the incenter of a triangle
    A=pts(1,:);
    B=pts(2,:);
    C=pts(3,:);
    D=pts(4,:); 
    
    incenter = (A + B + C + D)/4;
end