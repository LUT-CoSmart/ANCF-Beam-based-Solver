function [ux, uy, uz]= Displ(DofsAtNode,uu,nn)

    ux = uu(xlocANCF(DofsAtNode,nn,1)); 
    uy = uu(xlocANCF(DofsAtNode,nn,2));
    uz = uu(xlocANCF(DofsAtNode,nn,3));    
