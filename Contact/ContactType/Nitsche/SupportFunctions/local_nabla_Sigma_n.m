function dSigma_dxk = local_nabla_Sigma_n(Sigma_F,nabla_F,N)
    
    nabla_Sigma = Sigma_F * nabla_F; 
    dSigma_dxk = zeros(3,3); 
    
    for k = 1:3 
        dSigma_dxk(:,k) = reshape(nabla_Sigma(:,k),3,3) * N; 
    end
    