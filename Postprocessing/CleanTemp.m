function CleanTemp(Body, fl)

    if fl
        % delete body folder 
        BodyFolder = Body.BodyFolder;
        rmpath(BodyFolder);
        delete(fullfile(BodyFolder, '*'));
        rmdir(BodyFolder);
    
        % go to Temp
        folder = fileparts(BodyFolder);
    
        if isfolder(folder)
            contents = dir(folder);
            contents = contents(~ismember({contents.name}, {'.', '..'}));
    
            if isempty(contents)
                rmdir(folder);
            end
        end
    end