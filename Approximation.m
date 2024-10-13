% Gathering integration points
addpath('GaussPoints');
if App==0
   if (Area == 0) || (Area == 3) || (Area == 4)
      disp('****** The area type is not recognized ******');
      return;
   elseif Area == 1
      fac=1;
      ElementName = num2str(Element);             % using 'abcd' classification, see in https://doi.org/10.1007/s11071-022-07518-z
      num = str2double(ElementName(2));           % Number of point in axial direction, usually it equals to number of nodes
      [xi,wxi] = gauleg2(-1,1,num);
      pcirc(:,1)=repmat(xi',1,num);
      pcirc(:,2)=reshape(repmat(xi,1,num)',num^2,1);
      wvec=wxi';
      wcirc=reshape(wvec.*wvec',num^2,1);
   elseif Area == 2 
      fac=pi;
      pcirc=[[0,0]; [sqrt(2/3)*1,0]; [-sqrt(2/3)*1,0];[sqrt(1/6)*1,1/2*sqrt(2)];[sqrt(1/6)*1,-1/2*sqrt(2)];[-sqrt(1/6)*1,1/2*sqrt(2)];[-sqrt(1/6)*1,-1/2*sqrt(2)]];
      wcirc=[1/4,1/8,1/8,1/8,1/8,1/8,1/8];
   end   
else
   fac=1;
   [pcirc,wcirc]=PoiGen(data,nu2,App,Pointpic);   
end    
detF0=1/4*H*W*fac;
% Reorginizing points for AceGen
[Gint,Nint] = generateGint(n_xi,pcirc,wcirc); 
    