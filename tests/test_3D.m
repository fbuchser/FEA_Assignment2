%test_3D
nodes = [1 0.0 0.0 0.0; 2 10.0 0.0 0.0; 3 5.0 10*sin(60*pi/180) 0.0; 4 5.0 10/2*sin(60*pi/180) 10.0]
elements = [1 1 2 3 4 10000 0.2]

load1 = [4 100 200 100]
%load2 = []

bcs = [1 0 0 0; 2 0 0 0; 3 0 0 0]

[U, epsi, epsi1, epsi3, sigma, sigma1, sigma3] = solver_solid3D(nodes, elements, load1, bcs);
disp('Node 4 displacement:')
disp(U(10:12))