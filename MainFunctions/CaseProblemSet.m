function [Body,Force,Boundary] = CaseProblemSet(Body,CaseName,ApproximationScheme)
    
        
    switch CaseName
           case {"LockingBendingLarge", "LockingBendingSmall"}
                CSName = "Rectangular";
                Body.Length.X = 2; 
                Body.Length.Y = 0.5;
                Body.Length.Z = 0.1;
                MaterialName = 'KS';
                param.E = 2.07 * 1e11;
                param.nu = 0.3;                                 
                Force.Position.X = Body.Length.X;
                Boundary.Position = [];  
                Boundary.Type = "full";
                if CaseName == "LockingBendingSmall"
                   Force.Maginutude.Y = 62.5 * 1e3; 
                else
                   Force.Maginutude.Y = 62.5 * 1e6;
                end   

           case "ElongationHyperelastic"
                CSName = "Rectangular";
                Body.Length.X = 1; 
                Body.Length.Y = 0.1;
                Body.Length.Z = 0.1;
                MaterialName = 'Neo';                
                param.mu= 9 * 1e5;
                Force.Maginutude.X = 9500;
                Force.Position.X = Body.Length.X; 
                Boundary.Position = [];
                Boundary.Type = "reduced"; 
           
           case {'PrincetonBeamSmall', 'PrincetonBeamLarge'}
                CSName = "Rectangular";
                Body.Length.X = 0.508; 
                Body.Length.Y = 12.377 * 1e-3;
                Body.Length.Z = 3.2024 * 1e-3;
                MaterialName = 'KS';
                param.E = 71.7 * 1e9;
                param.nu = 0.31;                                 
                Force.Position.X = Body.Length.X;
                Boundary.Position = [];  
                Boundary.Type = "full";
                if CaseName == "PrincetonBeamSmall"
                   Force.Maginutude.Total = 8.896; 
                else 
                   Force.Maginutude.Total = 13.345;
                end 



           otherwise
                error('****** Case is not recognized ******');         
    end
    
    compressiblility= {'KS'};
    fibers= {'GOH'};

    Body = MaterialType(Body,MaterialName,param, compressiblility, fibers);
   
    addpath('GaussPoints');
    Body = GausPointsApprox(Body,CSName,ApproximationScheme);
   