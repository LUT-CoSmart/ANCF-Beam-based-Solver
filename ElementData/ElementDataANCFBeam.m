function Body = ElementDataANCFBeam(Body,ElmentType,SubType,Element) % Define element properties

Body.Element = Element;

if Element == 3243
    Body.Slope_x = true;  % Does the element include slope/derivative along x axis  
elseif (Element == 3333) || (Element==3343) || (Element==3353) || (Element==3363) || (Element==34103) 
    Body.Slope_x = false;
else    
    error('Element of type %s of subtype  %s is not recognized !!!!\n', ElmentType,SubType)
end 

ElementName = num2str(Element); % using 'abcd' classification, see in https://doi.org/10.1007/s11071-022-07518-z
Body.ElementNodes = str2double(ElementName(2));                % Number of nodes            
Body.DIM = str2double(ElementName(end));                       % Problem dimensionality     
Body.DofsAtNode = Body.DIM * str2double(ElementName(3:end-1)); % Number of Dofs in each element node
Body.ElementDofs=Body.ElementNodes*Body.DofsAtNode;            % Total number of Dofs in the chosen element
Body.ElementName = string(ElementName);

