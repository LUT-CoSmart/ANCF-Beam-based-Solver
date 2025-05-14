function Body = DefineElement(Body,ElementType,SubType,ElementName,Modification)
    
    addpath('ElementData')
    availbletypes = {'Beam'};

    if ismember(ElementType,availbletypes)
        switch SubType
               case "ANCF"  
                    if Modification == "None" 
                       Body = ElementDataANCFBeam(Body,ElementType,SubType,ElementName);
                    else
                       error('****** The modification is not available for the this Element Type (%s) and subtype (%s) ******',ElementType, SubType);  
                    end
             
        otherwise
               error('****** Element is not recognized of SubType %s of %s type  ******',SubType,ElementType);
        end
    else
        error('****** Element Type is not recognized ******');
    end
    
    Body.ElementType = ElementType;
    Body.SubType = SubType;
    
    

