% Subtendon Data for Surface 1
scale=1e-3;

curve1=scale*[
6.66 3.85;
6.41 4.10;
6.08 4.36;
5.55 4.74;
5.28 4.96;
5.05 5.21;
4.86 5.48;
4.67 5.84;
4.43 6.39;
4.42 6.41;
4.14 7.10;
];

curve2=scale*[
4.14 7.10;
4.01 7.50;
3.94 7.88;
3.89 8.34;
3.90 8.70;
3.94 8.98;
4.02 9.21;
];

curve3=scale*[
4.02 9.21;
4.13 9.42;
4.29 9.59;
4.48 9.74;
4.81 9.93;
5.02 10.01;
];

curve4=scale*[
5.02 10.01;
5.24 10.05;
5.44 10.04;
5.62 9.99;
5.77 9.90;
];

curve5=scale*[
5.77 9.90;
5.91 9.77;
6.00 9.67;
6.29 9.36;
6.48 9.21;
6.68 9.11;
6.91 9.05;
7.27 9.00;
];

curve6=scale*[
7.27 9.00;
7.57 8.95;
7.72 8.90;
7.89 8.79;
8.01 8.66;
];

curve7=scale*[
8.01 8.66;
8.13 8.47;
8.21 8.23;
8.26 7.96;
8.28 7.68;
8.29 7.37;
8.31 6.78;
8.38 5.54;
8.43 4.68;
8.43 4.21;
8.39 3.87;
];

curve8=scale*[
8.39 3.87;
8.35 3.72;
8.27 3.49;
8.18 3.35;
8.05 3.23;
];

curve9=scale*[
8.05 3.23;
7.87 3.16;
7.69 3.14;
7.47 3.16;
7.24 3.25;
];

curve10=scale*[
7.24 3.25;
7.04 3.38;
6.88 3.55;
6.66 3.85;
];

% collect all original data
%data_1=[curve1;curve2;curve3;curve4;curve5;curve6;curve7;curve8;curve9;curve10];
data_1={flipud(curve10);flipud(curve9);flipud(curve8);flipud(curve7);flipud(curve6);flipud(curve5);flipud(curve4);flipud(curve3);flipud(curve2);flipud(curve1)};

% max_x=max(data_1(:,1));
% min_x=min(data_1(:,1));
% max_y=max(data_1(:,2));
% min_y=min(data_1(:,2));
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
% 
% % collect changed data
% data_2=[curve1a;curve2a;curve3a;curve4a;curve5a;curve6a;curve7a;curve8a;curve9a;curve10a];
% 
% data={curve1a,curve2a,curve3a,curve4a,curve5a,curve6a,curve7a,curve8a,curve9a,curve10a};
% nu2=max(size(data));
