% Subtendon Data for Surface 1
scale=1e-3;

curve1=scale*[
6.48 14.01;
6.78 13.87;
7.11 13.68;
7.37 13.48;
7.61 13.25;
7.75 13.07;
7.87 12.87;
];

curve2=scale*[
7.87 12.87;
7.94 12.67;
8.00 12.41;
8.03 12.16;
8.12 11.09;
8.14 10.77;
8.19 10.20;
8.21 10.07;
8.25 9.76;
8.25 9.61;
8.21 9.47;
];

curve3=scale*[
8.21 9.47;
8.13 9.33;
7.99 9.18;
7.85 9.08;
7.66 9.00;
];

curve4=scale*[
7.66 9.00;
7.48 8.97;
7.28 8.98;
6.97 9.02;
6.73 9.08;
6.54 9.16;
6.36 9.29;
];

curve5=scale*[
6.36 9.29;
6.19 9.46;
5.96 9.72;
5.80 9.88;
5.65 9.97;
5.47 10.03;
];

curve6=scale*[
5.47 10.03;
5.30 10.04;
5.12 10.01;
4.89 9.93;
4.61 9.79;
4.43 9.73;
];

curve7=scale*[
4.43 9.73;
4.29 9.72;
4.14 9.75;
4.01 9.81;
];

curve8=scale*[
4.01 9.81;
3.92 9.91;
3.85 10.03;
3.80 10.18;
3.78 10.39;
];

curve9=scale*[
3.78 10.39;
3.79 10.70;
3.84 11.16;
3.90 11.49;
3.99 11.76;
4.12 12.00;
4.32 12.30;
4.48 12.58;
4.61 12.87;
];

curve10=scale*[
4.61 12.87;
4.76 13.31;
4.87 13.58;
4.98 13.75;
5.10 13.89;
5.19 13.97;
5.34 14.04;
];

curve11=scale*[
5.34 14.04;
5.50 14.09;
5.74 14.11;
6.03 14.12;
6.26 14.08;
6.48 14.01;
];

% collect all original data
%data_1={curve1;curve2;curve3;curve4;curve5;curve6;curve7;curve8;curve9;curve10;curve11};

data_1={flipud(curve11);flipud(curve10);flipud(curve9);flipud(curve8);flipud(curve7);flipud(curve6);flipud(curve5);flipud(curve4);flipud(curve3);flipud(curve2);flipud(curve1)};

% data_1s=[curve1;curve2;curve3;curve4;curve5;curve6;curve7;curve8;curve9;curve10;curve11];
% 
% max_x=max(data_1s(:,1));
% min_x=min(data_1s(:,1));
% max_y=max(data_1s(:,2));
% min_y=min(data_1s(:,2));
% 
% curve1a=[change(curve1(:,1),min_x,max_x) change(curve1(:,2),min_y,max_y)];
% curve2a=[change(curve2(:,1),min_x,max_x) change(curve2(:,2),min_y,max_y)];
% curve3a=[change(curve3(:,1),min_x,max_x) change(curve3(:,2),min_y,max_y)];
% curve4a=[change(curve4(:,1),min_x,max_x) change(curve4(:,2),min_y,max_y)];
% curve5a=[change(curve5(:,1),min_x,max_x) change(curve5(:,2),min_y,max_y)];
% curve6a=[change(curve6(:,1),min_x,max_x) change(curve6(:,2),min_y,max_y)];
% curve7a=[change(curve7(:,1),min_x,max_x) change(curve7(:,2),min_y,max_y)];
% curve8a=[change(curve8(:,1),min_x,max_x) change(curve8(:,2),min_y,max_y)];
% curve9a=[change(curve9(:,1),min_x,max_x) change(curve9(:,2),min_y,max_y)];
% curve10a=[change(curve10(:,1),min_x,max_x) change(curve10(:,2),min_y,max_y)];
% curve11a=[change(curve11(:,1),min_x,max_x) change(curve11(:,2),min_y,max_y)];
% 
% % collect changed data
% data_2={curve1a;curve2a;curve3a;curve4a;curve5a;curve6a;curve7a;curve8a;curve9a;curve10a;curve11a};
% 
% data={curve1a,curve2a,curve3a,curve4a,curve5a,curve6a,curve7a,curve8a,curve9a,curve10a,curve11a};
% nu2=max(size(data));
