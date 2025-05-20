function [data, nu2, CenterX, CenterY, CSSizeX, CSSizeY] = Binormalization(data_1)

    data_1_ver =  vertcat(data_1{:});

    maxPerColumn = max(data_1_ver);
    minPerColumn = min(data_1_ver);

    max_x = maxPerColumn(1);
    min_x = minPerColumn(1);
    max_y = maxPerColumn(2);
    min_y = minPerColumn(2);

    CenterX = (min_x+max_x)/2;
    CenterY = (min_y+max_y)/2;

    CSSizeX = max_x - min_x;
    CSSizeY = max_y - min_y;

    nu2=length(data_1);
    data = cell(nu2,1);
    for i = 1:nu2
        curve = data_1{i};
        data{i} = [change(curve(:,1),min_x,max_x) change(curve(:,2),min_y,max_y)];     
    end   

    Points = vertcat(data{:}); 

    % Point orientation checking
    x = Points(:,1);
    y = Points(:,2);
    direction  = 0.5 * sum(x(1:end-1) .* y(2:end) - x(2:end).*y(1:end-1));

    if direction < 0
       warning('PoiGen-used cross-sectional points are defined in clockwise order. Integration and Contact functions might give the opposite sign!');
    end 



