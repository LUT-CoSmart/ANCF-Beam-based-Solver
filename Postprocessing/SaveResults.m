function Body = SaveResults(Body, step, something)

    % Validate input
    if (nargin <3) ||...
        ~(isnumeric(something) || (isstring(something) && (something == "last" || something == "all"))) 
        
        warning("Provide correct identifier: must be numeric, 'last', or 'all', set to all");
        something = "all";
    end

    % Initialize or reset results
    if ~isfield(Body, 'Results') || (isstring(something) && something == "last")
        Body.Results = [];
    end

    % Default: do not save
    flag = false;

    % Save logic
    if isstring(something) && something == "all"
        flag = true;
    elseif isnumeric(something) && mod(step, something) == 0
        flag = true;
    elseif isstring(something) && something == "last"
        flag = true;  % this only works if called at the last step externally
    end

    % Save data
    if flag
        uf = Body.u(Body.fextInd);
        Body.Results = [Body.Results; Body.ElementNumber, Body.TotalDofs, uf'];
    end
