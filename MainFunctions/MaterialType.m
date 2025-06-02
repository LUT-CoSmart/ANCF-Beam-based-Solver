function  Body = MaterialType(Body,MaterialName,param)
    % Material parameters
    compressibleMaterials = {'KS'};
    fiberMaterials = {'GOH'};

    % Define "bulk" module 
    if ismember(MaterialName,compressibleMaterials)
       d = [];
    else
       d = 10^(-12); 
    end 

    % Define fibers
    if ismember(MaterialName,fiberMaterials)
       a0 = [1 0 0];   % fiber direction 
       Body.FiberTwist = 0; % inner (fiber) pre-twist
       Body.Fibers = true;
    else  % material isotrtopic
       a0 = [];
       Body.Fibers = false;
    end
    
    Body.MaterialName = MaterialName;
    Body.const = [cell2mat(struct2cell(param))', a0, d];


