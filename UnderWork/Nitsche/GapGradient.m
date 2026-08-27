function Outcome = GapGradient(Body, Outcome)

    n = numel(Outcome); % number of points of interest
    h = 1e-11;    
    
    % Contact knowledge
    Gaps = vertcat(Outcome.Gap);
    Normals = vertcat(Outcome.NormalProjected);    
    ContactPoints = vertcat(Outcome.ContactPoint);

    %% TODO: Make projection of all points at once, not 3 times 

    % rows = reshape(1:3*n, n, 3); % helping array to choose only one set of points for one xi coord 
    % ContactPoints_variation = repmat(ContactPoints, 3, 1); %  total array of points after variation   
    % 
    % for i = 1:3 % taking one full set (row(:,i) ) and i-th coord  
    %     ContactPoints_variation(rows(:,i),i) = ContactPoints_variation(rows(:,i),i) - h;
    % end
    % 
    % [Gaps_variation, ~, ~, Normals_variation, ~, ~] = ProjectPoints(Body, ContactPoints_variation);
    % 
    % Gaps_variation = reshape(Gaps_variation(:), n, 3);
    % Normals_variation = reshape(Normals_variation, n, 3, 3);


    for j = 1:n % Initialize new structure fields
        Outcome(j).nabla_Gap     = zeros(3,1);
        Outcome(j).nabla_Normals = zeros(3,3);
    end
    
    %% TODO: Maybe there is a sense to do via local coordinates 
    for i = 1:3
        
        ContactPoints_variation = ContactPoints; % Reset all coordinates 
        ContactPoints_variation(:,i) = ContactPoints(:,i) - h;

        % assumption that inside points stay inside after variation 
        [Gap_variation,~,~,Normals_variation,~,~]=ProjectPoints(Body,ContactPoints_variation);

        for j = 1:n
            Outcome(j).nabla_Gap(i) = (Gaps(j) - Gap_variation(j)) / h; 
            Outcome(j).nabla_Normals(i,:)= (Normals(j,:)' - Normals_variation(j,:)') / h;
        end    
        
    end 

    