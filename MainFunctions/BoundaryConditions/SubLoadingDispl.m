function Body = SubLoadingDispl(Body, currentStep, Nsteps, type)
    
    if currentStep == 1
       Prevloadstep = 0;
    else
       Prevloadstep = Body.Loadstep; 
    end

    switch type

           case "linear"
                Loadstep =  currentStep/Nsteps;
           case "exponential"     
                coef = 2;
                Loadstep = (exp((currentStep/Nsteps).^coef) - 1) / (exp(1) - 1);
           case "quadratic"
                Loadstep = (currentStep/Nsteps)^2;
           case "cubic"
                Loadstep = (currentStep/Nsteps)^3;        
           case "quartic"  
                Loadstep = (currentStep/Nsteps)^4;       
           case "logarithmic"                
                x = currentStep / Nsteps;
                Loadstep = log(1 + x) / log(2);
                
           case "mixed_Stepvise"
                threshold = 0.5; % threshold by number of steps
                power = 2;

                x = currentStep / Nsteps;   % normalized step [0,1]                                                              
                if x <= threshold  
                    Loadstep = x^power; % cubic part                  
                else
                    y_switch = threshold^power;   % value at switch
                    slope = (1 - y_switch) / (1 - threshold); % slope of linear part
                    Loadstep = y_switch + slope * (x - threshold);
                end   

           case "mixed_Loadvise"
                threshold = 0.5; % threshold by loads
                power = 3;

                x = currentStep / Nsteps;   % normalized step [0,1]                                                        
                x_switch = threshold^(1/power);   % x where cubic reaches threshold                
                if x <= x_switch
                    Loadstep = x^power;   % cubic part
                else
                    slope = (1 - threshold) / (1 - x_switch); % slope of linear part
                    Loadstep = threshold + slope * (x - x_switch);
                end
               
           otherwise
                error('Unknown loading type')                 
    end  

    SubDisplacement.Maginutude.X = Body.DisplacementVector(1) * (Loadstep-Prevloadstep);
    SubDisplacement.Maginutude.Y = Body.DisplacementVector(2) * (Loadstep-Prevloadstep);
    SubDisplacement.Maginutude.Z = Body.DisplacementVector(3) * (Loadstep-Prevloadstep);

    Body.applied_disp(Body.DisplInd) =  [SubDisplacement.Maginutude.X; SubDisplacement.Maginutude.Y; SubDisplacement.Maginutude.Z];
    Body.Loadstep = Loadstep;