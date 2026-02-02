function  Body = MaterialType(Body,MaterialName,param, compressiblility, fibers)
    
    Body.Volume =  Body.Length.X *  Body.Length.Y  *  Body.Length.Z;
    
    % Define "bulk" module 
    if nargin < 4 || ~ismember(MaterialName, compressiblility)
        %d = 1e-13;
        d = 1e-10 * Body.Volume; % emperical relation for incokmp. bodies with dependency of its size
    else
        d = [];
    end

    % Define fibers
    if nargin < 5 || ~ismember(MaterialName, fibers) % material isotrtopic
       param.a0 = [];
       Body.Fibers = false;
    else       
       Body.Fibers = true;
    end

    Body.MaterialName = MaterialName;
    Body.const = [cell2mat(struct2cell(param))', d];

   

    