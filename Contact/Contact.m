function [Kc,Fc,Gap,GapMax] = Contact(Body1,Body2,ContactType,ContactVariable,ContactRegType,ContactFiniteDiference)
    
    if ContactType == "None"
       Fc = zeros(Body1.TotalDofs + Body2.TotalDofs,1);
       Kc = zeros(length(Fc));
       Gap = NaN;
       GapMax.gap = 0;
       GapMax.area = NaN;
    else 
        %% TODO: add boxing to identify the necessity of the contact, for now we always consider its existence
        addpath("Contact\ContactType\");
        if ContactType == "Penalty"
           ContactType = @Penalty;
        elseif ContactType == "NitscheLin"
           ContactType = @NitscheLin; 
        else
           error('****** Contact type is not implemneted ******')
        end
        
        % Collecting bodies' surface points
        Body1.SurfacePoints = feval(Body1.SurfacefunctionName, Body1, Body1.q);         
        Body2.SurfacePoints = feval(Body2.SurfacefunctionName, Body2, Body2.q);
        
        switch ContactFiniteDiference
            case "Matlab"
                [Kc,Fc,Gap,GapMax] = ContactForceMatlab(Body1,Body2,ContactType,ContactVariable);

            case "Matlab_automatic"
                [Kc,Fc,Gap,GapMax] = ContactForceMatlabAuto(Body1,Body2,ContactType,ContactVariable);

            otherwise
                disp("Unknown finite difference scheme for the contact, switched to Matlab_automatic")
                [Kc,Fc,Gap,GapMax] = ContactForceMatlabAuto(Body1,Body2,ContactType,ContactVariable);
            
        end    
        [~,Kc] = Regularization(Kc,Fc,ContactRegType,false);
        
    end      
    
    