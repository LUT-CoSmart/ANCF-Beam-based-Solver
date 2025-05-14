function Body = Geometry(Body,CSName,ApproximationScheme)
    
        
    switch CSName
           case "Rectangular"
                Body.Length.X = 2; % Body length
                Body.Length.Y = 0.5;
                Body.Length.Z = 0.5;

           case "Oval"
                Body.Length.X = 2;   % Body length
                Body.Length.Y = 0.5; % Oval diameter in y-axis
                Body.Length.Z = 0.5; % Oval diameter in z-axis

           case {"C", "Tendon"}
                ApproximationScheme = "Poigen";
                disp("For chosen area the approximation scheme switched to Poigen")
                if CSName == "C"
                   Body.Length.X = 1; 
                   Body.Length.Y = 0.1; % Extension in y-axis
                   Body.Length.Z = 0.1; % Extension in z-axis
                elseif CSName == "Tendon"
                   Body.Length.X = 1;
                   Body.Length.Y = 0.1; % Extension in y-axis
                   Body.Length.Z = 0.1; % Extension in z-axis
                end                
                
           otherwise
                error('****** Cross-section is not recognized ******');
    end    
    
    addpath('GaussPoints');  
    Body = GausPointsApprox(Body,CSName,ApproximationScheme);
    

    

    
    