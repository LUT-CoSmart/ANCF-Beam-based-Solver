function Body = SurfaceBeamApprox(Body)    
    
    % Here are parameters to approximate Body surfaces in all directions
    Netazeta = 4; % used in standard CS (for Oval x4 multiplied)
    Nxi = 8;  % used in PoiGen & Standard

    ApproximationScheme = Body.IntegrationType;
    CSName = Body.CSName;
    IsoData = [];

    switch ApproximationScheme

           case "Poigen" % fixed CS 

               run(CSName); 
               [data, ~, ~, ~, ~, ~] = Binormalization(data_1);
               surfaceCSZetaEta = vertcat(data{:}); 

           case "Standard"                
               
                if CSName == "Rectangular"
                   % Surface data for contact and visualization 
                   % surfaceZetaEta = [-1  -1; 0  -1; 1  -1; 1   0; 1   1; 0   1; -1   1; -1   0; -1  -1];
                   % we need to create more dense discritization for the better projections
                   base = linspace(-1, 1, Netazeta);

                   bottom = [base', -1*ones(length(base),1)];
                   right  = [ones(length(base),1), base'];
                   top    = [flip(base)', ones(length(base),1)];
                   left   = [-1*ones(length(base),1), flip(base)'];
                   surfaceCSZetaEta = [bottom(1:end-1,:); right(1:end-1,:); top(1:end-1,:); left];
                   
                 elseif CSName == "Oval" 
                   % Surface data for contact and visualization
                   theta = linspace(0, 2*pi, 4*Netazeta);
                   surfaceCSZetaEta = [cos(theta)' sin(theta)'];
                 end
        
                 % Point orientation checking
                 eta = surfaceCSZetaEta(:,1);
                 zeta = surfaceCSZetaEta(:,2);
                 direction  = 0.5 * sum(eta(1:end-1) .* zeta(2:end) - eta(2:end).*zeta(1:end-1));
                
                 if direction < 0
                    warning('The surface points are defined in clockwise order. Contact functions might give the opposite sign!');
                 end  
                                  
        otherwise                  
                error('****** Unkown Integration type for %s ******', Body.Name);
                               
    end


    switch ApproximationScheme

           case {"Standard", "Poigen"} % fixed CS 
                surfaceCSZetaEta = unique(surfaceCSZetaEta, 'rows', 'stable');
                              
                pointCS = surfaceCSZetaEta;
                NpointCS = size(surfaceCSZetaEta,1);
                xi =  linspace(-1, 1, Nxi); 
                % Collect all points' isoparametric data 
                for k = 1:Body.ElementNumber                     
                    for i = 1:Nxi-1 % all but last layer per element (they are repeating in the next one)
                        for j = 1:NpointCS                                             
                            IsoData = [IsoData; xi(i), pointCS(j,2), pointCS(j,1), k];                             
                        end
                    end   
                end    
                % last layer to close the form
                for j = 1:NpointCS                
                    IsoData = [IsoData; xi(end), pointCS(j,2), pointCS(j,1), k];
                end

                 % Surface data for faces
                Body.SurfaceXi.Nxi = Nxi;
                Body.SurfaceXi.pointCS = surfaceCSZetaEta;

           otherwise                  
                error('****** Unkown Integration type for %s ******', Body.Name);  
    end

    % Surface data for contact and visualization       
    Body.IsoData = IsoData;
       
