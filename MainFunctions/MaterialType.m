function  Body = MaterialType(Body,MaterialName,param, compressiblility, fibers)

    % Define "bulk" module 
    if nargin < 4 || ~ismember(MaterialName, compressiblility)
        d = 1e-12;
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

   

    