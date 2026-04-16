function Sigma_F = SigmaJacobian(F_,fun)

    Sigma_vec = @(x) reshape( fun(reshape(x,3,3)),9, 1);
    y0 = F_(:);
    Sigma_vec0 = Sigma_vec(y0);
    G = @(t,y) Sigma_vec(y);
    fac = 1e-4;
    Sigma_F = numjac(G, 0, y0, Sigma_vec0, fac, []);

  