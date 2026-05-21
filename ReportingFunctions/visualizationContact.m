function visualizationContact(Gap,Body1,Body2,Show)
    

if Show  

    figure();
    visualization(Body1,Body1.q,'none',Show);    
    visualization(Body2,Body2.q,'none',Show);

    P = Gap.points;

    x = P(:,1);
    y = P(:,2);
    z = P(:,3);
    c = abs(P(:,4));

    tri = delaunay(z,y);  
    trisurf(tri, x, y, z, c, 'FaceColor','interp', 'EdgeColor','none');

    colormap jet
    colorbar
    title(['Total gap: ' ,num2str(Gap.total) ,'  , Maximum gap ', num2str(Gap.maximum)],'FontName','Times New Roman','FontSize',20);  

    plot3(x,y,z,'or') % checking contact points


else
    disp('No contact visualization')
end     