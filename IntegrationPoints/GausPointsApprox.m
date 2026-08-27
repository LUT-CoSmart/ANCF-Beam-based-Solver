function Body = GausPointsApprox(Body,CSName,ApproximationScheme, IntegrationPoints)    

    % Number of point in axial direction usually equals to number of nodes
    Body.IntegrationType = ApproximationScheme;    
    Body.CSName = CSName;
    num = Body.ElementNodes;
    
    if num < 3 % elemets with 2 nodes cannot recover stresses at all  
       num = 3; 
    end

    if IntegrationPoints == "Lobatto" % choosing initial for a main axis
        [xi,wxi] =  lglnodes(num-1);
    else    
        [xi,wxi] = gauleg2(-1,1,num);
    end  

    if ApproximationScheme == "Poigen"
       addpath("CrossSections") 
       run(CSName); 
       % Deg=input('Input Approximation degree (1 or above): ');      % Approximation degree for Green's formula    
       Deg = 1;    
       if  (CSName == "Rectangular") || (CSName == "RectangularCheck") ||(CSName== "Oval") % standard shapes
           [data, nu2, Body.CSCenterZ] = Binormalization(data_1); 
       else
           [data, nu2, Body.CSCenterZ, Body.CSCenterY, Body.Length.Z, Body.Length.Y] = Binormalization(data_1);              
       end
       [pcirc,wcirc]=PoiGen(data,nu2,Deg); 

    elseif ApproximationScheme == "Standard"

           if CSName == "Rectangular"  || (CSName == "RectangularCheck")                        
              pcirc(:,1)=repmat(xi',1,num);
              pcirc(:,2)=reshape(repmat(xi,1,num)',num^2,1);
              wvec=wxi';
              wcirc=reshape(wvec.*wvec',num^2,1);
                            
           elseif CSName == "Oval" 
              prefac =pi;              
              pcirc=[[0,0]; [sqrt(2/3)*1,0]; [-sqrt(2/3)*1,0];[sqrt(1/6)*1,1/2*sqrt(2)];[sqrt(1/6)*1,-1/2*sqrt(2)];[-sqrt(1/6)*1,1/2*sqrt(2)];[-sqrt(1/6)*1,-1/2*sqrt(2)]];
              wcirc=prefac*[1/4,1/8,1/8,1/8,1/8,1/8,1/8];       
           end
           
    elseif ApproximationScheme == "Volume"
         error('****** The approximation code was not added yet. ******');
    else
         error('****** Provide correct approximation code ******');
    end

    Body.detF0=1/4*Body.Length.Y*Body.Length.Z;
    
    % Reorginizing points for AceGen    
    Gint = []; % Initialize Gint matrix
    % Loop through all combinations of xiv, etav, and zetav
    for i = 1:length(xi)
        for j=1:length(pcirc)
            Gint= [Gint; xi(i), pcirc(j,1), pcirc(j,2), wxi(i) * wcirc(j)];
        end
    end
    Body.Nint = size(Gint,1); % total number of integration points 
    Body.Gint = Gint;
    Body.Area = (Body.Length.Y * Body.Length.Z / 8 ) *sum(Body.Gint(:,4));   
       
    