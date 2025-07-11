function Body = SubLoading(Body, currentStep, Nsteps, type, initial)
    
    if  nargin < 5 % preload
        initial = zeros(3,1);
    end 

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
    
    
    Subforce.Maginutude.X = (Body.ForceVectorInit(1) - initial(1)) * Loadstep;
    Subforce.Maginutude.Y = (Body.ForceVectorInit(2) - initial(2)) * Loadstep;
    Subforce.Maginutude.Z = (Body.ForceVectorInit(3) - initial(3)) * Loadstep;

    ForceVector = [Subforce.Maginutude.X; Subforce.Maginutude.Y; Subforce.Maginutude.Z];
    Fext(fextInd) = ForceVector;

    Body.Fext = Fext;