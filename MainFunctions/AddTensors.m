function Body = AddTensors(Body)
    
    path = 'TensorDerivations\' + Body.ElementType + '\' ;
    

    % Creation of temporal folder   
    TempRoot = fullfile(pwd, 'Temp');
    % Create Temp root if missing
    if ~isfolder(TempRoot)
        mkdir(TempRoot);
    end
    
    % Create folder for the specific body
    bodyFolder = fullfile(TempRoot, Body.Name);
    
    if ~isfolder(bodyFolder)
        mkdir(bodyFolder);
    end
    
    % Add shape function
    ShapeFunctionFolder = fullfile('TensorDerivations', Body.ElementType, 'Matlab', 'ShapeFunctions', Body.ElementName);
    copyfile(fullfile(ShapeFunctionFolder, '*'), bodyFolder, 'f');  % Copy all files from source to bodyFolder
    
    % Add files for inner energy calculations based on AceGen
    pathAceGen = path + 'AceGen';
    
    % Construct the function name
    function_name = 'ANCF'+Body.ElementName+Body.MaterialName;
            
    % Search recursively
    files = dir(fullfile(pathAceGen, '**', function_name + '.m') );
        
    % Check 
    if ~isempty(files)
        srcFile = fullfile(files(1).folder, files(1).name);
        destFile = fullfile(bodyFolder, 'AceGenForce.m');
        copyfile(srcFile, destFile, 'f'); 
            
        % some modifications to AceGen file for MEX (surprisingly it also speed up AceGen)
        lines = readlines(destFile);
            
        pattern_check = contains(lines(15), "persistent v") && contains(lines(16), "if size(v)<") && contains(lines(17), "v=zeros");
            
        if pattern_check
           % Extract the number N from 'v=zeros(N,...'
           n = regexp(lines(17), 'v\s*=\s*zeros\((\d+)', 'tokens', 'once');
           if ~isempty(n)
               lines(15) = "v=zeros(" + n{1} + ",'double');";
               lines(16) = "";
               lines(17) = "";
               lines(18) = "";
           end
         end
         writelines(lines, destFile);
           
         else
           warning('This element is not yet implemented in AceGen, substitude by the dummy');
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           fileName = 'AceGenForce.m';
           fullPath = fullfile(bodyFolder, fileName);
           dummyCode = [
                    "function [q0f, q0, uu, D, K, F, GINT, nintpt] = AceGenForce(q0f, q0, uu, D, K, F, GINT, nintpt)"
                    "K = zeros(" + Body.ElementDofs + ");"
                    "F = zeros(" + Body.ElementDofs + ",1);"
                    "end"];
            fid = fopen(fullPath, 'w');
            for i = 1:numel(dummyCode)
                fprintf(fid, '%s\n', dummyCode(i));
            end
            fclose(fid);     
     end

     % Add files for inner energy calculations based on Matlab
     srcFolder = fullfile(pwd, 'MaterialsSecondKirhhoff');
     srcFile = fullfile(srcFolder, Body.MaterialName + ".m"); 
     destFile = fullfile(bodyFolder, 'PiolaSecondTensor.m'); % Material check has already done in 'Materials.m'
     copyfile(srcFile, destFile, 'f'); % force overwrite

     pathMatlab = path + 'Matlab' + '\' + Body.SolutionBase;
     if Body.SolutionBase == "Position"
        pathMatlab =  pathMatlab + '\' + Body.ElementName;
     else
        pathMatlab =  pathMatlab +  '\' + Body.DeformationType + '\' + Body.ElementName;
     end
     
     files = dir(fullfile(pathMatlab, '*.*'));
     files = files(~[files.isdir]);
     
     for k = 1:length(files)
        srcFile = fullfile(pathMatlab, files(k).name);
        destFile = fullfile(bodyFolder, files(k).name);
        copyfile(srcFile, destFile, 'f');  % 'f' to force overwrite
     end


    switch Body.FiniteDiference 
       case "AceGen"   
            disp("For chosen finite difference scheme, deformations are ony finite and displacement-based") 
    
            Body.DeformationType = "Finite";
            Body.SolutionBase = "Displacement";
                    
       case "Matlab"   

            if Body.SolutionBase == "Position"
               disp("For chosen finite difference and solution-based scheme, deformations are only finite")
               Body.DeformationType = "Finite"; 
            end
    
       otherwise
            error('****** Choose correct Finite Diference scheme ******\n')
    end        

  
    
    addpath(bodyFolder);
    Body.BodyFolder = bodyFolder;