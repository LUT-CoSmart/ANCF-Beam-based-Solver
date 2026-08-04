function Body = AddTensors(Body)
        
    path = 'TensorDerivations\' + Body.ElementType + '\' ;
    
    % Creation of temporal folder   
    
    TempRoot = fullfile(pwd, 'Temp');
    % Create Temp root if missing
    if ~isfolder(TempRoot)
        mkdir(TempRoot);
    else 

    end
    
    % Create folder for the specific body
    bodyFolder = fullfile(TempRoot, Body.Name);
    
    if ~isfolder(bodyFolder)
        mkdir(bodyFolder);
    end
    
    % Add shape function
    ShapeFunctionFolder = fullfile('TensorDerivations', Body.ElementType, 'Matlab', 'ShapeFunctions', Body.ElementName);
    copyfile(fullfile(ShapeFunctionFolder, '*'), bodyFolder, 'f');  % Copy all files from source to bodyFolder
    
    % 
    pathAceGen = path + 'AceGen';    
    pathMatlab = path + 'Matlab\ElementFunctions' + '\' + Body.ElementName;

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
           warning('This element is not yet implemented in AceGen, substituded by the dummy');
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
          
     files = dir(fullfile(pathMatlab, '*.*'));
     files = files(~[files.isdir]);
     
     for k = 1:length(files)
        srcFile = fullfile(pathMatlab, files(k).name);
        destFile = fullfile(bodyFolder, files(k).name);
        copyfile(srcFile, destFile, 'f');  % 'f' to force overwrite
     end

     % To be sure, the folders with files of the same names are removed  
     rmpath(genpath(fullfile(pwd, 'TensorDerivations')));

     switch Body.FiniteDiference 
       case "AceGen"   
            disp("For chosen finite difference scheme, deformations are ony finite and displacement-based")    
            Body.DeformationType = "Finite";
            Body.SolutionBase = "Displacement";
                    
       case {"Matlab", "Matlab_automatic"}   

            if Body.SolutionBase == "Position"
               disp("For chosen finite difference and solution-based scheme, deformations are only finite")
               Body.DeformationType = "Finite"; 
            end
           
       otherwise
            error('****** Choose correct Finite Diference scheme ******\n')
     end        
    
     addpath(bodyFolder);    
     Body.BodyFolder = bodyFolder;

     Body.SurfacefunctionName = "Build" + Body.ElementType + "Surface"; 
        
     L = Body.Length.Ln;
     W = Body.Length.Z;
     H = Body.Length.Y;
    
     Body.Shape = @(xi,eta,zeta) Shape_(L,H,W,xi,eta,zeta);       
     
     Shape_xi = @(xi,eta,zeta) Shape_xi_(L,H,W,xi,eta,zeta);
     Shape_eta = @(xi,eta,zeta) Shape_eta_(L,H,W,xi,eta,zeta);
     Shape_zeta = @(xi,eta,zeta) Shape_zeta_(L,H,W,xi,eta,zeta);

     Body.Nm_xi_ = @(xi,eta,zeta) [Shape_xi(xi,eta,zeta); Shape_eta(xi,eta,zeta); Shape_zeta(xi,eta,zeta)];
     Body.nabla_r_xi = @(xi,eta,zeta,q) [Shape_xi(xi,eta,zeta)*q Shape_eta(xi,eta,zeta)*q Shape_zeta(xi,eta,zeta)*q];        
     
     Body.F0 = @(q0,xi,eta,zeta) F0(q0,L,H,W,xi,eta,zeta);
     Body.F = @(q,q0,u,xi,eta,zeta) DeformationGradient(Body.SolutionBase,q,q0,u,L,H,W,xi,eta,zeta);   

     if Body.Fibers
        a0 = Body.Dvec(end-6:end-4)';
        Body.a0_fib = @(q0,Phi,xi,eta,zeta) a0_fib(a0,q0,Phi,L,H,W,xi,eta,zeta)';
     end
        
     if Body.Fibers
        S = @(F_, a0_fib) PiolaSecondTensor(F_, Body.const, a0_fib);
     else
        S = @(F_) PiolaSecondTensor(F_, Body.const);         
     end

     Sigma = @(F_,S) 1/det(F_) * F_ * S * F_'; %  Cauchy Stresses 
    
     Body.S = S;
     Body.Sigma = Sigma;
     Body.Sigma_nn = @(F_, N, S) N' * Sigma(F_, S) * N;    

           
     Body.NodeSphere = feval("MaxNode" + Body.ElementType + "Dimension", Body); % space around node for possible contact check;
     
     % if Body.mex 
     %    CreateMex(create,Body);
     %    InnerForce = @(Body) InnerForce_mex(Body);
     % end

     Body.Results = [];


     % For contact: in this way of the points' organization, the normals of the trimesh is directed to the inside volume 
     % Option: 
     % SurfacePoints = BuildBeamSurface(Body,Body.q0);
     % faces = Body.BodyFaces;
     % [mean_nodes,face_normals]=getFaceCenterAndNormals(faces,SurfacePoints);
     % quiver3(mean_nodes(:,1), mean_nodes(:,2), mean_nodes(:,3),face_normals(:,1),  face_normals(:,2),  face_normals(:,3), 0.5, 'r', 'LineWidth', 0.01); 
     % patch('Vertices',SurfacePoints,'Faces',faces,'FaceColor','cyan','EdgeColor','black');