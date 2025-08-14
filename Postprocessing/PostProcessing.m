function PostProcessing(Body,Results,visDeformed,visInitial)
    disp('***********************************************')
    fprintf('Static test for %s, Element %s \n', Body.ElementType, Body.ElementName)
    fprintf('Deformations = %s \n', Body.DeformationType)
    fprintf('n & DOFs & ux & uy & uz \n')
    for k=1:size(Results,1) 
      fprintf('%d & %d & %10.8f & %10.8f & %10.8f  \n',Results(k,1:5))
    end

    fprintf('Reactions = %s \n', Body.DeformationType)
    fprintf('Rx & Ry & Rz \n')
    for k=1:size(Results,1) 
      fprintf('%10.8f & %10.8f & %10.8f  \n',Results(k,6:8))
    end
    
    
    if (visDeformed~=true) && (visInitial~=true) 
        disp('No visualization')
        return
    else
        figure;
        % axis equal
        set(gca, 'FontSize', [12], 'FontName','Times New Roman');
        set(text, 'FontSize', [12], 'FontName','Times New Roman');
        xlabel('\it{X}','FontName','Times New Roman','FontSize',[20])
        ylabel('\it{Y}','FontName','Times New Roman','FontSize',[20]),
        zlabel('Z [m]','FontName','Times New Roman','FontSize',[20]);
        grid minor
        
        visualization(Body,Body.q0,'red',visInitial);
        visualization(Body,Body.q,'cyan',visDeformed);
    end    

