function [Kc,Fc,Gap] = Contact(Body1,Body2,ContactTypeName,ContactVariable,ContactFiniteDiference,CalculateStiffness)
    
    if nargin < 6 % some methods don't need to calculate stiffness matrix
        CalculateStiffness = true;
    end
    
    Fc = zeros(Body1.TotalDofs + Body2.TotalDofs,1);
    Kc = zeros(length(Fc));
    Gap.total = NaN;
    Gap.maximum= 0;
    Gap.area = NaN;     
    Gap.points = [];

    if ContactTypeName ~= "None" 
        %% TODO: add boxing to identify the necessity of the contact, for now we always consider its existence
        switch ContactTypeName
            
            case "Penalty"
                ContactType = @Penalty;
        
            case "Nitsche"
                ContactType = @Nitsche;
        
            otherwise
                error("Contact type '%s' is not implemented.", ContactTypeName);
        end
        
        % Collecting bodies' surface points        
        Body1.SurfacePoints = Body1.SurfacePointsFunction(Body1.q);       
        Body2.SurfacePoints = Body2.SurfacePointsFunction(Body2.q);

        switch ContactFiniteDiference
            case "Matlab"
                [Kc,Fc,Gap] = ContactForceMatlab(Body1,Body2,ContactType,ContactVariable,CalculateStiffness);

            case "Matlab_automatic"
                [Kc,Fc,Gap] = ContactForceMatlabAuto(Body1,Body2,ContactType,ContactVariable,CalculateStiffness);

            otherwise
                disp("Unknown finite difference scheme for the contact, switched to Matlab_automatic")
                [Kc,Fc,Gap] = ContactForceMatlabAuto(Body1,Body2,ContactType,ContactVariable,CalculateStiffness);          
        end    
    end      
    
    