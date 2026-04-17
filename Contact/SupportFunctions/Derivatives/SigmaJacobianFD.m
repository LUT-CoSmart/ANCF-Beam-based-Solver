function Sigma_F = SigmaJacobianFD(F_, fun)
    Sigma_F = zeros(9,9);
    h = 1e-6;

    for i = 1:3
        for j = 1:3
            F_plus = F_;
            F_minus = F_;

            F_plus(i,j) = F_plus(i,j) + h;
            F_minus(i,j) = F_minus(i,j) - h;

            Sigma_plus = fun(F_plus);
            Sigma_minus = fun(F_minus);

            col = sub2ind([3,3], i, j);
            Sigma_F(:,col) = (Sigma_plus(:) - Sigma_minus(:)) / (2*h);
        end
    end
end