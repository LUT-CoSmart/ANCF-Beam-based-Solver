function Body = Geometry(Body,CSName,ApproximationScheme)
    
        
    switch CSName
           case "Rectangular"
                Body.Length.X = 2; % Body length
                Body.Length.Y = 0.5;
                Body.Length.Z = 0.5;
                Body.Volume =  Body.Length.X *  Body.Length.Y  *  Body.Length.Z;
                
           case "Oval"
                Body.Length.X = 38e-3;   % Body length
                Body.Length.Y = 0.00655; % Oval diameter in y-axis
                Body.Length.Z = 0.00655; % Oval diameter in z-axis

           case {"C", "Tendon",...
                 "Middle_cross_section1_1", "Middle_cross_section2_1", "Middle_cross_section3_1",...
                 "Sol_subj2_middle", "MG_subj2_middle", "LG_subj2_middle",...
                 "ten_Sol_3", "ten_MG_3", "ten_LG_3",...
                 "ten_Sol_3_2", "ten_MG_3_2", "ten_LG_3_2"}
                ApproximationScheme = "Poigen";
                disp("For chosen area the approximation scheme switched to Poigen")
                
                if (CSName == "C") || (CSName == "Tendon")
                   Body.Length.X = 1;
                elseif (CSName == "Middle_cross_section1_1") || (CSName == "Middle_cross_section2_1") || (CSName == "Middle_cross_section3_1") || ...
                        (CSName == "ten_Sol_3") || (CSName == "ten_MG_3") || (CSName == "ten_LG_3")  || ...
                        (CSName == "ten_Sol_3_2") || (CSName == "ten_MG_3_2") || (CSName == "ten_LG_3_2")  
                   Body.Length.X = 0.07; 
                elseif (CSName == "Sol_subj2_middle") || (CSName == "MG_subj2_middle") || (CSName == "LG_subj2_middle") 
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

    


    
    

    
    