%% Subtendon
scale=1/1000*0.0536;      
curve1=[17  30
81  30
137 30]*scale;
curve2=[137 30
140 33
142 35
144 38
145 41
140 48
131 57
125 64
120 71
119 75
118 76
120 80
125 87
133 95
135 98
137 99]*scale;
curve3=[137 99
146 95
160 90
171 83
181 75
187 65
190 58
191 52
190 47
186 35
176 25
159 14
142  7
122  3
108  2
104  1
99   0
94   1
83   2
73   3
63   5
52   9
46  10
37  14
31  18
23  23
17  30]*scale;
% collect all original data
data_1=[curve1;curve2;curve3]; 

max_x=max(data_1(:,1));
min_x=min(data_1(:,1));
max_y=max(data_1(:,2));
min_y=min(data_1(:,2));

%% Such as all point in clockwise direction (see lines 84-90); we need to reorder them
curve1a=[change(curve1(:,1),min_x,max_x) change(curve1(:,2),min_y,max_y)];
curve1a=flipud(curve1a);
curve2a=[change(curve2(:,1),min_x,max_x) change(curve2(:,2),min_y,max_y)];
curve2a=flipud(curve2a);
curve3a=[change(curve3(:,1),min_x,max_x) change(curve3(:,2),min_y,max_y)];
curve3a=flipud(curve3a);
data_2=[curve3a;curve2a;curve1a];
data={curve3a,curve2a,curve1a};
nu2=max(size(data));