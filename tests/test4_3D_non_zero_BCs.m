%test4_3D_non_zero_BCs
nodes = [1 0.0 0.0 0.0;
    2 10.0 0.0 0.0;
    3 10.0 10.0 0.0;
    4 0.0 10.0 0.0;
    5 0.0 0.0 10.0;
    6 10.0 0.0 10.0;
    7 10.0 10.0 10.0;
    8 0.0 10.0 10.0
    ];
elements = [1 1 2 3 6 10000 0.2;
    2 1 3 4 8 10000 0.2;
    3 5 8 6 1 10000 0.2;
    4 6 7 3 8 10000 0.2];

%load1 = [5 0 0 -25;
%    6 0 0 -25;
%    7 0 0 -25;
%    8 0 0 -25];
load1 = [6 0 0 -25;
    7 0 0 -25];
%load2 = []

%load1 = [5 -10 -10]
%load2 = []

bcs = [1 0 0 0; 2 0 0 0; 3 0 0 0; 4 0 0 0; 5 0 0 0; 6 0.001 0.001 NaN; 7 0.001 0.001 NaN; 8 0 0 0]

[U, epsi, epsi1, epsi3, sigma, sigma1, sigma3] = solver_solid3D(nodes, elements, load1, bcs);
disp('Node 6 displacement (should be ~0.001, 0.001, free):')
disp(U(16:18))
disp('Node 7 displacement (should be ~0.001, 0.001, free):')
disp(U(19:21))