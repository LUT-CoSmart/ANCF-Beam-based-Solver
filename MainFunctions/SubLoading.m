function Body = SubLoading(Body, currentStep, Nsteps, type)
    
    Fext = Body.Fext;
    fextInd = Body.fextInd;
    
    switch type
           case "linear"
                Loadstep =  currentStep/Nsteps;
           case "exponential"     
                coef = 5;
                Loadstep = (exp((currentStep/Nsteps).^coef) - 1) / (exp(1) - 1);
           case "quadratic"
                Loadstep = (currentStep/Nsteps)^2;
           case "cubic"
                Loadstep = (currentStep/Nsteps)^3;    
           case "quartic"  
                Loadstep = (currentStep/Nsteps)^4;    
          
           otherwise
                error('Unknown loading type')                 
    end  
    
    
    Subforce.Maginutude.X = Body.ForceVectorInit(1) * Loadstep;
    Subforce.Maginutude.Y = Body.ForceVectorInit(2) * Loadstep;
    Subforce.Maginutude.Z = Body.ForceVectorInit(3) * Loadstep;

    ForceVector = [Subforce.Maginutude.X; Subforce.Maginutude.Y; Subforce.Maginutude.Z];
    Fext(fextInd) = ForceVector;

    Body.Fext = Fext;