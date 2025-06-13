function Body = Geometry(Body,CSName,ApproximationScheme,k)
    
        
    switch CSName
           case "Rectangular"
                % Body.Length.X = 1; % Body length
                % Body.Length.Y = 0.1;
                % Body.Length.Z = 0.1;

                Body.Length.X = 43e-3*k; 
                Body.Length.Y = 2.14028e-05*k;
                %Body.Length.Y = 0.01;
                Body.Length.Z = 1*k;
                
                Body.Volume =  Body.Length.X *  Body.Length.Y  *  Body.Length.Z;
           case "Oval"
                Body.Length.X = 2;   % Body length
                Body.Length.Y = 0.5; % Oval diameter in y-axis
                Body.Length.Z = 0.5; % Oval diameter in z-axis

           case {"C", "Tendon",...
                 "Middle_cross_section1_1", "Middle_cross_section2_1", "Middle_cross_section3_1",...
                 "Sol_subj2_middle"}
                ApproximationScheme = "Poigen";
                disp("For chosen area the approximation scheme switched to Poigen")
                
                if (CSName == "C") || (CSName == "Tendon")
                   Body.Length.X = 1;
                elseif (CSName == "Middle_cross_section1_1") || (CSName == "Middle_cross_section2_1") || (CSName == "Middle_cross_section3_1")   
                   Body.Length.X = 0.07; 
                elseif (CSName == "Sol_subj2_middle")   
                   Body.Length.X = 43e-3; 
                end                
                
        otherwise

               filepath = fullfile("CrossSections", CSName + ".m"); 
               if exist(filepath, 'file') == 2
                  ApproximationScheme = "Poigen";
                  disp("For chosen area the approximation scheme switched to Poigen") 
               else
                  error('****** Cross-section is not recognized ******');
               end  
    end    
    
    addpath('GaussPoints');  
    Body = GausPointsApprox(Body,CSName,ApproximationScheme);

    
    

    
    