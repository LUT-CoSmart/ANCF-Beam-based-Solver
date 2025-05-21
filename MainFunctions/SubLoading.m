function Subforce = SubLoading(Force,currentStep, Nsteps, type)

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
    
    % Force check
    if ~isfield(Force.Maginutude, 'X')
       Force.Maginutude.X = 0;
    end
    if ~isfield(Force.Maginutude, 'Y')
       Force.Maginutude.Y = 0;
    end
    if ~isfield(Force.Maginutude, 'Z')
       Force.Maginutude.Z = 0;
    end 

    % Force position check
    if ~isfield(Force.Position, 'X')
       Force.Position.X = 0;
    end
    if ~isfield(Force.Position, 'Y')
       Force.Position.Y = 0;
    end
    if ~isfield(Force.Position, 'Z')
       Force.Position.Z = 0;
    end 


    Subforce.Maginutude.X = Force.Maginutude.X * Loadstep;
    Subforce.Maginutude.Y = Force.Maginutude.Y * Loadstep;
    Subforce.Maginutude.Z = Force.Maginutude.Z * Loadstep;
    Subforce.Position = Force.Position;