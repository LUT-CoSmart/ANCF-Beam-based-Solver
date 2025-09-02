%% Subtendon (LG, Type 3)
scale=1/3000*0.0536;

curve1=scale*[
        330 167
        329 260
        327 352
        % 330 167
        % 317 969
        % 308 1524
        % 316 1030
        % 327 352
        ];

curve2=scale*[327 352
        233 347
        128 317
        40 257
        8 187
        13 156
        26 130];
 
curve3 = scale * [  
    % 330 167
    % 187 150
    % 26 130
    330 167
    102 140
    26 130
];

curve3=flipud(curve3);

data_1={curve1;curve2;curve3};
