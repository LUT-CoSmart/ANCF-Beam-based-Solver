
function CreateMex(create,Body)

writeMex = false;
if create
    origFolder = pwd;
    cd InnerForceFunctions;
    
    mexFilePath = fullfile('InnerForceFunctions', 'InnerForce_mex');
    if exist(mexFilePath, 'file') == 3
        disp('InnerForce_mex already exists.');
        answer = input('Do you want to rewrite it? (y/n): ', 's');
        if lower(answer) == 'y'           
            writeMex = true;
        end    
    else 
        writeMex = true;
    end
    
    if writeMex
            %% TODO: check does these changes affect the body behaviour
            % they connect to folders , where numerical data are comming
            % from
            
            % Body.ElementName = coder.typeof('a', [1, Inf], [false, true]);
            % Body.ElementType = coder.typeof('a', [1, Inf], [false, true]);
            % Body.SubType = coder.typeof('a', [1, Inf], [false, true]);
            % Body.DeformationType = coder.typeof('a', [1, Inf], [false, true]);
            % Body.FiniteDiference = coder.typeof('a', [1, Inf], [false, true]);
            % Body.SolutionBase = coder.typeof('a', [1, Inf], [false, true]);
            % Body.Name = coder.typeof('a', [1, Inf], [false, true]);            
            % Body.IntegrationType = coder.typeof('a', [1, Inf], [false, true]);
            % Body.CSName = coder.typeof('a', [1, Inf], [false, true]);
            % Body.BodyFolder = coder.typeof('a', [1, Inf], [false, true]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Body.MaterialName = coder.typeof('a', [1, Inf], [false, true]);
            Body.const = coder.typeof(0, [1, Inf], [false, true]);  
            Body.Gint = coder.typeof(0, [Inf, 4], [true, false]);  
            Body.const = coder.typeof(0, [1, Inf], [false, true]);  
            Body.PosDofs = coder.typeof(0, [1, Inf], [false, true]);
            Body.nloc = coder.typeof(0, [Inf, Inf], [true, true]);
            Body.phim = coder.typeof(0, [Inf, Inf], [true, true]);
            Body.Phim = coder.typeof(0, [Inf, Inf], [true, true]);
            Body.xloc = coder.typeof(0, [Inf, Inf], [true, true]);
            Body.u = coder.typeof(0, [Inf, 1], [true, false]);
            Body.q0 = coder.typeof(0, [Inf, 1], [true, false]);
            Body.q = coder.typeof(0, [Inf, 1], [true, false]);
            Body.q0f = coder.typeof(0, [Inf, 1], [true, false]);
            Body.RotInv = coder.typeof(0, [Inf, Inf], [true, true]);
            Body.Rot = coder.typeof(0, [Inf, Inf], [true, true]);
            Body.P00 = coder.typeof(0, [Inf, 3], [true, false]);
            Body.Dvec = coder.typeof(0, [1, Inf], [false, true]);
            Body.SurfaceXi.pointCS = coder.typeof(0, [Inf, 2], [true, false]);
            Body.IsoData = coder.typeof(0, [Inf, 4], [true, false]);
            Body.BodyFaces = coder.typeof(0, [Inf, 3], [true, false]);
            % 
            Body.bc = coder.typeof(false, [1, Inf], [false, true]);
            Body.Fext = coder.typeof(0, [Inf, 1], [true, false]);
            Body.fextInd = coder.typeof(0, [1, Inf], [false, true]);
            Body.ForceVectorInit = coder.typeof(0, [Inf, 1], [false, true]);

            codegen InnerForce -args {Body} -config:mex
    end
    
    if isfolder('codegen')
       rmdir('codegen', 's');
    end
    
    cd(origFolder);
end
