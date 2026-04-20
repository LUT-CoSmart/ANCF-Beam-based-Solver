function [Kc,Fc,Gap,GapMax] = Contact(Body1,Body2,ContactTypeName,ContactVariable,ContactFiniteDiference,CalculateStiffness)
    
    if nargin < 6 % some methods don't need to calculate stiffness matrix
        CalculateStiffness = true;
    end
    
    Fc = zeros(Body1.TotalDofs + Body2.TotalDofs,1);
    Kc = zeros(length(Fc));
    Gap = NaN;
    GapMax.gap = 0;
    GapMax.area = NaN;     
    
    if ContactTypeName ~= "None" 
        %% TODO: add boxing to identify the necessity of the contact, for now we always consider its existence
        % Penalty approach        
        if ContactTypeName == "Penalty"
           ContactType = @Penalty;

        % Nitshes approach           
        elseif contains(ContactTypeName, "Nitsche")
             ContactType = @(varargin)Nitsche(ContactTypeName,varargin{:});

        else
           error('****** Contact type is not implemneted ******')
        end
        
        % Collecting bodies' surface points
        Body1.SurfacePoints = feval(Body1.SurfacefunctionName, Body1, Body1.q);         
        Body2.SurfacePoints = feval(Body2.SurfacefunctionName, Body2, Body2.q);
        
        switch ContactFiniteDiference
            case "Matlab"
                [Kc,Fc,Gap,GapMax] = ContactForceMatlab(Body1,Body2,ContactType,ContactVariable,CalculateStiffness);

            case "Matlab_automatic"
                [Kc,Fc,Gap,GapMax] = ContactForceMatlabAuto(Body1,Body2,ContactType,ContactVariable,CalculateStiffness);

            otherwise
                disp("Unknown finite difference scheme for the contact, switched to Matlab_automatic")
                [Kc,Fc,Gap,GapMax] = ContactForceMatlabAuto(Body1,Body2,ContactType,ContactVariable,CalculateStiffness);          
        end    
    end      
    
    