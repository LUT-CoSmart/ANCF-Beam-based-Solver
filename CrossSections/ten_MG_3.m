%% Subtendon (MG, Type 3)
scale=1/3000*0.0536;
curve1 = scale * [
    607 130
    550 140
    474 150
    380 160
    330 167        
];


curve2 = scale * [
    330 167
    187 157
    26 130     
];

curve3 = scale * [
    26 130
    153 42
    232 22 
    394 20
    488 43          
    607 130];

%curve3 = scale * [
    % 26 130
    % 76 80
    % 153 42
    % 232 22
    % 316 14   
    % 394 20
    % 488 43
    % 563 80            
    % 607 130];


% collect all original data
data_1={curve1;curve2;curve3};
