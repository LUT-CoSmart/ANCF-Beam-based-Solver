function Body = GausPointsApprox(Body,CSName,ApproximationScheme)    

    % Number of point in axial direction usually equals to number of nodes
    Body.IntegrationType = ApproximationScheme;    
    Body.CSName = CSName;
    num = Body.ElementNodes;
        
    if ApproximationScheme == "Poigen"

       addpath("CrossSections") 
       run(CSName); 
       Deg=input('Input Approximation degree (1 or above): ');      % Approximation degree for Green's formula        
       [data, nu2] = Adjustment(data_1);
       [pcirc,wcirc]=PoiGen(data,nu2,Deg);
       prefac = 1;

    elseif ApproximationScheme == "Standard"

           if CSName == "Rectangular"
              prefac = 1;     
              [xi,wxi] = gauleg2(-1,1,num);
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

    Body.detF0=1/4*Body.Length.Y*Body.Length.Z*prefac;

    % Reorginizing points for AceGen
    [Body.Gint,Body.Nint] = generateGint(num,pcirc,wcirc); 
    
       
