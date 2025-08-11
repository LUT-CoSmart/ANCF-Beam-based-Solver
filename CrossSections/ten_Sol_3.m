%% Subtendon (Sol, Type 3)
scale=1/3000*0.0536;

% curve1=scale*[330 167
%         317 212
%         308 262
%         316 311
%         327 352];

curve1=scale*[330 167
        308 262
        327 352];

curve2=scale*[607 130
        550 140
        474 150
        380 160
        330 167
        ];
curve3=scale*[327 352
        444 339
        529 308
        594 260
        626 184
        623 158
        607 130];    

% curve3 = scale * [
%     607 130
%     400 160
%     380 165
%     330 167
% ];


% collect all original data

curve1=flipud(curve1);

curve2=flipud(curve2);

curve3=flipud(curve3);

data_1={curve1;curve2;curve3};
