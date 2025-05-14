function Body = Materials(Body, MaterialName)
            
        
    switch MaterialName
           case "Neo" % Neo-Hookean
                param.mu=9e5;       

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
                param.c10 = 7.64e3;
                param.k1 = 996.6e3;
                param.k2 = 524.6;      
                param.kappa = 0;   % fiber dipersion

           case "KS" % Kirhhoff-Saint-Venant
                param.E=2.07e11;
                param.nu=0.3;
                
           otherwise   
                error('****** The material type is not recognized ******');
    end
    

    Body = MaterialType(Body,MaterialName,param);
    
