function Body = Materials(Body, MaterialName, Subtype)
            
    if  nargin < 3 % preload
        Subtype = "";
    end 

    switch MaterialName
           case "Neo" % Neo-Hookean
                            
               switch  Subtype
                   case ""
                       param.mu=9e5;
                   
                   case "optimized_SEE"  
                       param.mu=28559355.9219739;

                   case "Recalculated_to_245mm"
                       param.mu=90e6/(0.245/0.07);    

                   case "Alex"
                       param.mu=90e6;
                    
                   case "Sol_old"
                       param.mu=103.1e6;

                   case "MG_old"
                       param.mu=143.2e6;

                   case "LG_old"
                       param.mu=226.7e6; 
               end

           case "Mooney2" % 2 contant Mooney-Rivlin
                param.c10=33.4e4;
                param.c01 = -337;

           case "Mooney5" % 5-contant Mooney-Rivlin 
                param.c10 = -7.7e5;
                param.c01 = 9.1e5;
                param.c11 = 1.03e6;
                param.c20 = -2.7e5;
                param.c02 = -5.9e5;
                 
           case "GOH" % Gasser-Ogden-Holzaphel material  

               switch  Subtype
                   case ""
                       param.c10 = 7.64e3;
                       param.k1 = 996.6e3;
                       param.k2 = 524.6;      
                       param.kappa = 0;      % fiber dipersion
                       param.a0 = [1 0 0]';   % fiber direction 
                       Body.FiberTwist = 0; % inner (fiber) pre-twist
    
                   case "Alex"
                       param.c10 = 53600;
                       param.k1 = 7.5351e7;
                       param.k2 = 23.926;      
                       param.kappa = 0;      % fiber dipersion
                       param.a0 = [1 0 0]';   % fiber direction 
                       Body.FiberTwist = 0; % inner (fiber) pre-twist
                   
                   case "Amir"
                        param.c10 = 9.67e6;
                        param.k1 = 135.5e6;
                        param.k2 = 131;
                        param.kappa = 0; % fiber dipersion
                        param.a0 = [1 0 0]'; % fiber direction
                        Body.FiberTwist = 0; % inner (fiber) pre-twist
                         
                   case "optimized_SEE"    
                       param.c10 = 53600;
                       param.k1 = 8409291.05908193;
                       param.k2 = 4.78018491589778;      
                       param.kappa = 0;      % fiber dipersion
                       param.a0 = [1 0 0]';   % fiber direction 
                       Body.FiberTwist = 0; % inner (fiber) pre-twist 
                   
                   case "stiff"
                       param.c10 = 50e6;
                       param.k1 = 7.5351e7;
                       param.k2 = 23.926;      
                       param.kappa = 0;      % fiber dipersion
                       param.a0 = [1 0 0]';   % fiber direction 
                       Body.FiberTwist = 0; % inner (fiber) pre-twist 

               end
                 
   
           case "KS" % Kirhhoff-Saint-Venant
                param.E=2.07e11;
                param.nu=0.3;
                
           otherwise   
                error('****** The material type is not recognized ******');
    end
    

    compressiblility= {'KS'};
    fibers= {'GOH'};

    Body = MaterialType(Body,MaterialName, param, compressiblility, fibers);
    
    