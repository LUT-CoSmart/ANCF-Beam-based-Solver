function  Body = AssigningMaterialParameters(Body,MaterialName,param, compressiblility, fibers)
    
    Body.Volume =  Body.Length.X *  Body.Length.Y  *  Body.Length.Z;
    
    % Define "bulk" module 
    if ~ismember(MaterialName, compressiblility)
        % d = 1e-13;
        d = 1e-8 * Body.Volume; % emperical relation for incokmp. bodies with dependency of its size
    else
        d = [];
    end

    % Define fibers
    if ~ismember(MaterialName, fibers) % material isotrtopic
       Body.Fibers = false;
    else       
       Body.Fibers = true;
    end

    Body.MaterialName = MaterialName;
    Body.const = [cell2mat(struct2cell(param))', d];

   

    