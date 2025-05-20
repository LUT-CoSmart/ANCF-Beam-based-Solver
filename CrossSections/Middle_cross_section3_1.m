% Subtendon Data for Surface 1
scale=1;

curve1=[
5.30 2.40;
5.48 2.00;
5.73 1.46;
5.84 1.27;
5.99 1.08;
6.29 0.78;
6.55 0.56;
6.80 0.38;
7.07 0.25;
7.37 0.13;
];

curve2=[
7.37 0.13;
7.64 0.05;
7.94 0.01;
8.35 0.01;
8.71 0.03;
8.90 0.06;
9.07 0.11;
9.20 0.19;
];

curve3=[
9.20 0.19;
9.30 0.28;
9.38 0.40;
9.43 0.57;
9.46 0.80;
];

curve4=[
9.46 0.80;
9.47 0.87;
9.45 1.11;
9.38 1.46;
9.31 1.71;
9.22 1.92;
9.10 2.12;
8.85 2.49;
8.72 2.70;
8.62 2.85;
];

curve5=[
8.62 2.85;
8.51 2.95;
8.37 3.03;
8.18 3.08;
7.86 3.13;
7.51 3.18;
7.29 3.24;
];

curve6=[
7.29 3.24;
7.14 3.30;
7.00 3.41;
6.88 3.54;
6.65 3.86;
6.44 4.09;
6.24 4.27;
];

curve7=[
6.24 4.27;
6.06 4.39;
5.87 4.48;
5.67 4.53;
5.49 4.54;
5.34 4.52;
];

curve8=[
5.34 4.52;
5.20 4.46;
5.09 4.39;
5.02 4.30;
];

curve9=[
5.02 4.30;
4.97 4.19;
4.96 4.04;
4.98 3.73;
5.00 3.56;
5.06 3.16;
5.16 2.79;
5.30 2.40;
];

% collect all original data
data_1=[curve1;curve2;curve3;curve4;curve5;curve6;curve7;curve8;curve9];

max_x=max(data_1(:,1));
min_x=min(data_1(:,1));
max_y=max(data_1(:,2));
min_y=min(data_1(:,2));

curve1a=[change(curve1(:,1),min_x,max_x) change(curve1(:,2),min_y,max_y)];
curve2a=[change(curve2(:,1),min_x,max_x) change(curve2(:,2),min_y,max_y)];
curve3a=[change(curve3(:,1),min_x,max_x) change(curve3(:,2),min_y,max_y)];
curve4a=[change(curve4(:,1),min_x,max_x) change(curve4(:,2),min_y,max_y)];
curve5a=[change(curve5(:,1),min_x,max_x) change(curve5(:,2),min_y,max_y)];
curve6a=[change(curve6(:,1),min_x,max_x) change(curve6(:,2),min_y,max_y)];
curve7a=[change(curve7(:,1),min_x,max_x) change(curve7(:,2),min_y,max_y)];
curve8a=[change(curve8(:,1),min_x,max_x) change(curve8(:,2),min_y,max_y)];
curve9a=[change(curve9(:,1),min_x,max_x) change(curve9(:,2),min_y,max_y)];

% collect changed data
data_2=[curve1a;curve2a;curve3a;curve4a;curve5a;curve6a;curve7a;curve8a;curve9a];

data={curve1a,curve2a,curve3a,curve4a,curve5a,curve6a,curve7a,curve8a,curve9a};
nu2=max(size(data));
