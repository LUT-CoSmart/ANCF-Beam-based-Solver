



origFolder = pwd;
cd InnerForceFunctions;

mexFilePath = fullfile('InnerForceFunctions', 'InnerForce_mex');
if exist(mexFilePath, 'file') == 3
    disp('InnerForce_mex already exists.');
    answer = input('Do you want to rewrite it? (y/n): ', 's');
    if lower(answer) == 'y'
        codegen InnerForce -args {Body} -config:mex
    end    
else
    codegen InnerForce -args {Body} -config:mex
    
end


if isfolder('codegen')
   rmdir('codegen', 's');
end

cd(origFolder);
