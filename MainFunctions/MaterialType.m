function  Body = MaterialType(Body,MaterialName,param)
    % Material parameters
    compressibleMaterials = {'KS'};
    fiberMaterials = {'GOH'};

    % Define "bulk" module 
    if ismember(MaterialName,compressibleMaterials)
       Body.d = [];
    else
       Body.d = 10^(-14); 
    end 

    % Define fibers
    if ismember(MaterialName,fiberMaterials)
       Body.a0 = [1 0 0];   % fiber direction 
       Body.FiberTwist = 0; % inner (fiber) pre-twist
       Body.Fibers = true;
    else  % material isotrtopic
       Body.a0 = [];
       Body.Fibers = false;
    end
    
    Body.MaterialName = MaterialName;
    Body.const = [cell2mat(struct2cell(param))', Body.a0, Body.d];


