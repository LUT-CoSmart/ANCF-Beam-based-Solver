function visualization_StressRecovery(Body,Show,name,FromGausElement)
    
    
    if Show 
        figure();
        axis equal
        visualization(Body,Body.q,'none',Show); 
        NameGaus = '';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        [Recovery_S1,Recovery_S2] = GeneralizeCoorStress(Body);

        xloc = Body.xloc;
        faces = Body.BodyFaces;
        IsoData = Body.IsoData;
        Shape   = Body.Shape;
        SurfacePoints = Body.SurfacePointsFunction(Body.q); 
           
        getFaceCenterAndNormals = Body.getFaceCenterAndNormals;
        [~,face_normals]=getFaceCenterAndNormals(faces,SurfacePoints);
        
        N_points = size(IsoData,1);
        points = [SurfacePoints zeros(N_points,1)];
        counter = zeros(N_points,1);
        
        if name == "Normal"
            Interested_stress = @(Sigma,Normal) Stress_nn(Sigma,Normal);
        elseif name == "VM"
            Interested_stress = @(Sigma,Normal) Stress_VM(Sigma);
        else
            error("This stress type is not implemneted!")
        end 
         
        if FromGausElement
            IsoData = GausElement(Body,IsoData); % recovering from Gaus element  
            NameGaus = ' from Gaus Element';
        end

        for k = 1:size(faces,1)  
            
            normal = face_normals(k,:)';
            Xi_points = IsoData(faces(k,:),:);
            n = size(Xi_points, 1); % Determine the number of surface points (3 or 4)
            Element = Xi_points(1,4); % all face points belong to one element        
            
            Generalize_S1_k = Recovery_S1(xloc(Element,:));
            Generalize_S2_k = Recovery_S2(xloc(Element,:));
                        
            for ii = 1:n          

                Stress_S1 = Shape(Xi_points(ii,1),Xi_points(ii,2),Xi_points(ii,3)) * Generalize_S1_k;
                Stress_S2 = Shape(Xi_points(ii,1),Xi_points(ii,2),Xi_points(ii,3)) * Generalize_S2_k;

                s_recovered = [Stress_S1; Stress_S2];

                Sigma = [s_recovered(1), s_recovered(4), s_recovered(6);
                         s_recovered(4), s_recovered(2), s_recovered(5);
                         s_recovered(6), s_recovered(5), s_recovered(3)];

                Stress_value  = Interested_stress(Sigma,normal);
                    
                point_number = faces(k,ii);
                points(point_number,4) = points(point_number,4) + Stress_value ;
                counter(point_number) = counter(point_number) + 1;
            end

        end
        points(:,4) = points(:,4)./counter;
        h = patch('Faces',faces,'Vertices',points(:,1:3),'FaceVertexCData',points(:,4),'FaceColor','interp','EdgeColor','k');

        colormap(turbo);
        cb = colorbar;
        cb.Label.String = [name ' stress '];
        title([name ' stress via Generalize stress coordinates and shape funciton' NameGaus]);

    else
        disp('No visualization')
    end     
end

function Stress = Stress_VM(Sigma)
        Sigma_dev = Sigma - trace(Sigma)/3 * eye(3); % Deviatoric stress tensor
        Stress = sqrt(3/2 * sum(Sigma_dev.^2, 'all')); % von Mises stress   
end
    
function Stress = Stress_nn(Sigma,Normal)
        Stress = Normal' * Sigma * Normal;
end
