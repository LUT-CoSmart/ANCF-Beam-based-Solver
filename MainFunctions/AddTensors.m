function Body = AddTensors(Body)
        
    path = 'TensorDerivations\' + Body.ElementType + '\' ;
    
    
    
    TempRoot = fullfile(pwd, 'Temp'); % Creation of temporal folder   
    % Create Temp root if missing
    if ~isfolder(TempRoot)
        mkdir(TempRoot);
    else 
        bodyFolder = fullfile(TempRoot, Body.Name);
        if isfolder(bodyFolder)
           rmdir(bodyFolder,'s');
        end   
    end
    
    % Create folder for the specific body
    bodyFolder = fullfile(TempRoot, Body.Name);
    
    if ~isfolder(bodyFolder)
        mkdir(bodyFolder);
    end
    
    % Add shape function
    ShapeFunctionFolder = fullfile('TensorDerivations', Body.ElementType, 'Matlab', 'ShapeFunctions', Body.ElementName);
    copyfile(fullfile(ShapeFunctionFolder, '*'), bodyFolder, 'f');  % Copy all files from source to bodyFolder

    % Construct the function name
    function_name = 'ANCF'+Body.ElementName+Body.MaterialName;
    
     switch Body.FiniteDiference 
       case "AceGen"   
            disp("For chosen finite difference scheme, deformations are ony finite and displacement-based")    
            Body.DeformationType = "Finite";
            Body.SolutionBase = "Displacement";
               

            pathAceGen = path + 'AceGen'; 
            files = dir(fullfile(pathAceGen, '**', function_name + '.m') );
                
            if ~isempty(files) % Check 
                srcFile = fullfile(files(1).folder, files(1).name);
                destFile = fullfile(bodyFolder, 'AceGenForce.m');
                copyfile(srcFile, destFile, 'f'); 
                    
                % some modifications to AceGen file for MEX
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
                error('This element is not yet implemented in AceGen');  
            end

       case {"Matlab", "Matlab_automatic"}   
             
                  
             pathMatlab = path + 'Matlab\ElementFunctions' + '\' + Body.ElementName;
             files = dir(fullfile(pathMatlab, '*.*'));
             files = files(~[files.isdir]);
             
             for k = 1:length(files)
                srcFile = fullfile(pathMatlab, files(k).name);
                destFile = fullfile(bodyFolder, files(k).name);
                copyfile(srcFile, destFile, 'f');  % 'f' to force overwrite
             end
            
             % Add files for inner energy calculations based on Matlab
             srcFolder = fullfile(pwd, 'MaterialsSecondKirhhoff');
             srcFile = fullfile(srcFolder, Body.MaterialName + ".m"); 
             destFile = fullfile(bodyFolder, 'PiolaSecondTensor.m'); % Material check has already done in 'Materials.m'
             copyfile(srcFile, destFile, 'f'); % force overwrite
             % To be sure, the folders with files of the same names are removed  
             


            if Body.SolutionBase == "Position"
               disp("For chosen finite difference and solution-based scheme, deformations are only finite")
               Body.DeformationType = "Finite"; 
            end
           
         case "Casadi"
            import casadi.*
            disp("For chosen finite difference scheme, deformations are ony finite and position-based")    
            Body.DeformationType = "Finite";
            Body.SolutionBase = "Position";
            pathCasadi= path + 'Casadi';
            files = dir(fullfile(pathCasadi, '**', function_name + '.casadi') );
            if ~isempty(files) % Check 
                srcFile = fullfile(files(1).folder, files(1).name);
                destFile = fullfile(bodyFolder, 'Casadi.casadi');
                copyfile(srcFile, destFile, 'f'); 
                Body.CasadiForce = Function.load(char(destFile));

            else
                error('This element is not yet implemented with Casadi');  
            end

             
       otherwise
            error('****** Choose correct Finite Diference scheme ******\n')
     end    

     rmpath(genpath(fullfile(pwd, 'TensorDerivations')));

     addpath(bodyFolder);    
     Body.BodyFolder = bodyFolder;

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Body.SurfacefunctionName = "Build" + Body.ElementType + "Surface";         
     L = Body.Length.Ln;
     W = Body.Length.Z;
     H = Body.Length.Y;
    
     Body.Shape = @(xi,eta,zeta) Shape_(L,H,W,xi,eta,zeta);       
     
     Shape_xi = @(xi,eta,zeta) Shape_xi_(L,H,W,xi,eta,zeta);
     Shape_eta = @(xi,eta,zeta) Shape_eta_(L,H,W,xi,eta,zeta);
     Shape_zeta = @(xi,eta,zeta) Shape_zeta_(L,H,W,xi,eta,zeta);

     Nm_xi_ = @(xi,eta,zeta) [Shape_xi(xi,eta,zeta); Shape_eta(xi,eta,zeta); Shape_zeta(xi,eta,zeta)];
     
     Body.nabla_r_xi = @(q,xi,eta,zeta) [Shape_xi(xi,eta,zeta)*q Shape_eta(xi,eta,zeta)*q Shape_zeta(xi,eta,zeta)*q];        
     Body.F = @(q,q0,u,xi,eta,zeta) DeformationGradient(Body.SolutionBase,q,q0,u,L,H,W,xi,eta,zeta);
     Body.F0 = @(q0,u,xi,eta,zeta) F0(q0,L,H,W,xi,eta,zeta);
     Body.dF_dq = @(q0,xi,eta,zeta) pagemtimes(reshape(Nm_xi_(xi,eta,zeta),3,3,[]), F0(q0,L,H,W,xi,eta,zeta)^(-1));
        
     % adjustment for vectorization
     if Body.DeformationType == "Small"
        Body.dEdq = @(dF,F) 0.5 * (pagetranspose(dF) + dF);
     else
        Body.dEdq = @(dF,F) 0.5 * (pagemtimes(dF,'transpose',F,'none') + pagemtimes(F,'transpose',dF,'none'));
     end

     if Body.Fibers 
        a0 = Body.Dvec(end-6:end-4)';
        a0_fiber = @(q0,Phi,xi,eta,zeta) a0_fib(a0,q0,Phi,L,H,W,xi,eta,zeta)';
        S = @(F_,q0,Phi,xi,eta,zeta) PiolaSecondTensor(F_, Body.const, a0_fiber(q0,Phi,xi,eta,zeta));
     else % in the nonfiber version, the arguments are accepted but simply ignored
        S = @(F_,q0,Phi,xi,eta,zeta) PiolaSecondTensor(F_, Body.const);         
     end

     Body.S = S; % 2nd Piola Stress tensor 
     Body.Sigma = @(F_,S) 1/det(F_) * F_ * S * F_'; %  Cauchy Stress tensor
           
     Body.NodeSphere = feval("MaxNode" + Body.ElementType + "Dimension", Body); % space around node for possible contact check;
     
     % if Body.mex 
     %    CreateMex(create,Body);
     %    InnerForce = @(Body) InnerForce_mex(Body);
     % end

     Body.Results = [];
     SurfaceShape = BuildBeamSurfaceMatrixLocal(Body);    
     Body.SurfacePointsFunction = @(q) reshape(SurfaceShape*q,3,[])';   
