function [u,epsi,epsi1,epsi3,sigma,sigma1,sigma3] = solver_solid3D(nodes,elements,load1,bcs,plot_scale)
%NUMBER of NODES
N = length(nodes(:,1));

%X- and y- coordinates of nodes
x = nodes(:,2);
y = nodes(:,3);
z = nodes(:,4);

%Number of elements
n = length(elements(:,1));

%Global stiffness matrix
K = zeros(3*N,3*N);
k1 = 0;
save('stiffness_matrices','k1')
for i = 1:n;
    E = elements(i,6);
    nyy = elements(i,7);

    
    %Nodes
    I = elements(i,2);
    J = elements(i,3);
    Q = elements(i,4);
    L = elements(i,5);
    
    x1 = nodes(I,2);
    x2 = nodes(J,2);
    x3 = nodes(Q,2);
    x4 = nodes(L,2);
    y1 = nodes(I,3);
    y2 = nodes(J,3);
    y3 = nodes(Q,3);
    y4 = nodes(L,3);
    z1 = nodes(I,4);
    z2 = nodes(J,4);
    z3 = nodes(Q,4);
    z4 = nodes(L,4);

    
    % Volume
    x21=x2-x1; x31=x3-x1; x41=x4-x1;
    y21=y2-y1; y31=y3-y1; y41=y4-y1;
    z21=z2-z1; z31=z3-z1; z41=z4-z1;
    V = (x21*(y31*z41 - y41*z31) + y21*(x41*z31 - x31*z41) + z21*(x31*y41 - x41*y31)) / 6;

    %Difference Variables & Coefficients
    y43=y4-y3; y42=y4-y2; y32=y3-y2; y41=y4-y1; y31=y3-y1; y21=y2-y1;
    z43=z4-z3; z42=z4-z2; z32=z3-z2; z41=z4-z1; z31=z3-z1; z21=z2-z1;
    x41=x4-x1; x31=x3-x1; x21=x2-x1;
    
    a1 =  y2*z43 - y3*z42 + y4*z32;
    a2 = -y1*z43 + y3*z41 - y4*z31;
    a3 =  y1*z42 - y2*z41 + y4*z21;
    a4 = -y1*z32 + y2*z31 - y3*z21;

    b1 = -x2*z43 + x3*z42 - x4*z32;
    b2 =  x1*z43 - x3*z41 + x4*z31;
    b3 = -x1*z42 + x2*z41 - x4*z21;
    b4 =  x1*z32 - x2*z31 + x3*z21;

    c1 =  x2*y43 - x3*y42 + x4*y32;
    c2 = -x1*y43 + x3*y41 - x4*y31;
    c3 =  x1*y42 - x2*y41 + x4*y21;
    c4 = -x1*y32 + x2*y31 - x3*y21;

    % B-Matrix
    BT = 1/(6*V) * [
        a1  0   0   a2  0   0   a3  0   0   a4  0   0 ;
        0  b1  0    0  b2  0    0  b3  0    0  b4  0 ;
        0   0  c1   0   0  c2   0   0  c3   0   0  c4;
        b1  a1  0   b2  a2  0   b3  a3  0   b4  a4  0 ;
        0  c1  b1   0  c2  b2   0  c3  b3   0  c4  b4;
        c1  0   a1  c2  0   a2  c3  0   a3  c4  0   a4];
    
    B = BT';
    % C-Matrix
    
    C = E/((1+nyy)*(1-2*nyy)) * [ 1-nyy  nyy    nyy    0              0              0;
                                  nyy    1-nyy  nyy    0              0              0;
                                  nyy    nyy    1-nyy  0              0              0;
                                  0      0      0      (1-2*nyy)/2    0              0;
                                  0      0      0      0              (1-2*nyy)/2    0;
                                  0      0      0      0              0              (1-2*nyy)/2];
    %Stiffness Matrix k
    k = V*B*C*B';
    
    if i == 1
        k1 = k;
        save('stiffness_matrices','k1','-append')
    end
    if i == 2
        k2 = k;
        save('stiffness_matrices','k2','-append')
    end
    if i == 3
        k3 = k;
        save('stiffness_matrices','k3','-append')
    end
    if i == 4
        k4 = k;
        save('stiffness_matrices','k4','-append')
    end
    
    
    %DOFs
    I1=3*I-2; I2=3*I-1; I3=3*I;
    J1=3*J-2; J2=3*J-1; J3=3*J;
    Q1=3*Q-2; Q2=3*Q-1; Q3=3*Q;
    L1=3*L-2; L2=3*L-1; L3=3*L;  
    index = [I1 I2 I3 J1 J2 J3 Q1 Q2 Q3 L1 L2 L3];
    
    %Global stiffness matrix
    K_fake = zeros(3*N,3*N);
    K_fake(index,index) = k;
    K = K+K_fake;
end


%Load vector
NN = size(load1,1);
rp = zeros(3*N,1);
for i = 1:NN,          %for each element
    f = load1(i,2:4); %load Fx,Fy,Fz
    I_x = 3*load1(i,1) - 2; % get X-DOF
    J_y = 3*load1(i,1) - 1; % get Y-DOF
    L_z = 3*load1(i,1) - 0; % get Z-DOF
    index = [I_x J_y L_z];
    rp(index) = f;
end
F = rp;

%Load - distributed 
%The load is specified in terms of element edge numbers. They are numbered
%1, 2 and 3
 % XXXX_------------------------ removed distributed loads
 % ------------------------

%BCs
NN = size(bcs,1);
F_old = F;
for i = 1:NN,
    r = bcs(i,2:4);
    I = 3*bcs(i,1) - 2;
    J = 3*bcs(i,1) - 1;
    L = 3*bcs(i,1) - 0;
    index = [I J L];
    F(index) = r;

    %The effect of BCs on the global stiffness matrix
    if ~isnan(r(1)), K(I,:)=0; K(I,I)=1; end
    if ~isnan(r(2)), K(J,:)=0; K(J,J)=1; end
    if ~isnan(r(3)), K(L,:)=0; K(L,L)=1; end
end

I_nan = isnan(F);
F(I_nan) = F_old(I_nan);
%pause

%Solution
u = K\F;


%Stress and strain
epsi = zeros(n,7);
epsi(:,1) = 1:n;
sigma = zeros(n,7);
sigma(:,1) = 1:n;

%Copy from above ot extract C and B
for i = 1:n,
    E = elements(i,6);
    nyy = elements(i,7);
    
    %Nodes
    I = elements(i,2);
    J = elements(i,3);
    Q = elements(i,4);    L = elements(i,5);
    
    I1 = 3*I - 2; I2 = 3*I - 1; I3 = 3*I - 0;
    J1 = 3*J - 2; J2 = 3*J - 1; J3 = 3*J - 0;
    Q1 = 3*Q - 2; Q2 = 3*Q - 1; Q3 = 3*Q - 0;
    L1 = 3*L - 2; L2 = 3*L - 1; L3 = 3*L - 0;
    index = [I1 I2 I3 J1 J2 J3 Q1 Q2 Q3 L1 L2 L3];
    d = u(index);
    
    %Nodes
    I = elements(i,2);
    J = elements(i,3);
    Q = elements(i,4);
    L = elements(i,5);
    
    x1 = nodes(I,2);
    x2 = nodes(J,2);
    x3 = nodes(Q,2);
    x4 = nodes(L,2);
    y1 = nodes(I,3);
    y2 = nodes(J,3);
    y3 = nodes(Q,3);
    y4 = nodes(L,3);
    z1 = nodes(I,4);
    z2 = nodes(J,4);
    z3 = nodes(Q,4);
    z4 = nodes(L,4);

    % Volume
    x21=x2-x1; x31=x3-x1; x41=x4-x1;
    y21=y2-y1; y31=y3-y1; y41=y4-y1;
    z21=z2-z1; z31=z3-z1; z41=z4-z1;
    V = (x21*(y31*z41 - y41*z31) + y21*(x41*z31 - x31*z41) + z21*(x31*y41 - x41*y31)) / 6;

    %Difference Variables & Coefficients
    y43=y4-y3; y42=y4-y2; y32=y3-y2; y41=y4-y1; y31=y3-y1; y21=y2-y1;
    z43=z4-z3; z42=z4-z2; z32=z3-z2; z41=z4-z1; z31=z3-z1; z21=z2-z1;
    x41=x4-x1; x31=x3-x1; x21=x2-x1;
    
    a1 =  y2*z43 - y3*z42 + y4*z32;
    a2 = -y1*z43 + y3*z41 - y4*z31;
    a3 =  y1*z42 - y2*z41 + y4*z21;
    a4 = -y1*z32 + y2*z31 - y3*z21;

    b1 = -x2*z43 + x3*z42 - x4*z32;
    b2 =  x1*z43 - x3*z41 + x4*z31;
    b3 = -x1*z42 + x2*z41 - x4*z21;
    b4 =  x1*z32 - x2*z31 + x3*z21;

    c1 =  x2*y43 - x3*y42 + x4*y32;
    c2 = -x1*y43 + x3*y41 - x4*y31;
    c3 =  x1*y42 - x2*y41 + x4*y21;
    c4 = -x1*y32 + x2*y31 - x3*y21;

    % B-Matrix
    BT = 1/(6*V) * [
        a1  0   0   a2  0   0   a3  0   0   a4  0   0 ;
        0  b1  0    0  b2  0    0  b3  0    0  b4  0 ;
        0   0  c1   0   0  c2   0   0  c3   0   0  c4;
        b1  a1  0   b2  a2  0   b3  a3  0   b4  a4  0 ;
        0  c1  b1   0  c2  b2   0  c3  b3   0  c4  b4;
        c1  0   a1  c2  0   a2  c3  0   a3  c4  0   a4];

    C = E/((1+nyy)*(1-2*nyy)) * [ 1-nyy  nyy    nyy    0              0              0;
        nyy    1-nyy  nyy    0              0              0;
        nyy    nyy    1-nyy  0              0              0;
        0      0      0      (1-2*nyy)/2    0              0;
        0      0      0      0              (1-2*nyy)/2    0;
        0      0      0      0              0              (1-2*nyy)/2];
    
    B = BT'; %12x6 for  k = V*B*C*B'
    eps = BT*d;

    epsi(i,2:7) = eps.';
    sigma(i,2:7) = (C*eps).';

    %PRINCIPLE STRESSES
    ex  = epsi(i,2); ey  = epsi(i,3); ez  = epsi(i,4);
    gxy = epsi(i,5); gyz = epsi(i,6); gzx = epsi(i,7);

    %Strain Tensor
    E_tensor = [ex      gxy/2   gzx/2;
        gxy/2   ey      gyz/2;
        gzx/2   gyz/2   ez   ];

    ev = sort(eig(E_tensor));
    epsi3(i) = ev(1);   % min principal strain
    epsi1(i) = ev(3);   % max principal strain

    % Stress tensor
    sx  = sigma(i,2); sy  = sigma(i,3); sz  = sigma(i,4);
    txy = sigma(i,5); tyz = sigma(i,6); tzx = sigma(i,7);

    S_tensor = [sx  txy tzx;
        txy sy  tyz;
        tzx tyz sz];

    sv = sort(eig(S_tensor));
    sigma3(i) = sv(1);  % min principal stress
    sigma1(i) = sv(3);  % max principal stress
    
end

%DRAW DRAW  (3D version)
n = length(elements(:,1));

% Coordinates
x = nodes(:,2);
y = nodes(:,3);
z = nodes(:,4);

min_x = min(x);
max_x = max(x);
min_y = min(y);
max_y = max(y);
min_z = min(z);
max_z = max(z);

b = 0.1 * max([mean([min_x max_x]), mean([min_y max_y]), mean([min_z max_z])]);

%% 1) DRAW INPUT: undeformed mesh
figure(1); clf;
tri = elements(:,2:5);          % 4-node tets
coords = [x y z];               % N x 3

tetramesh(tri, coords, 'FaceAlpha', 0.1, 'EdgeColor', 'k');
axis equal;
xlim([min_x-b, max_x+b]);
ylim([min_y-b, max_y+b]);
zlim([min_z-b, max_z+b]);
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
title('Undeformed mesh');

% OPTIONAL: label elements at centroids (can comment out if too busy)
NN = size(elements,1);
hold on;
for i = 1:NN
    I = elements(i,2);
    J = elements(i,3);
    Q = elements(i,4);
    L = elements(i,5);

    cx = mean(x([I J Q L]));
    cy = mean(y([I J Q L]));
    cz = mean(z([I J Q L]));

    h = text(cx, cy, cz, num2str(i), 'HorizontalAlignment','center', ...
                                   'VerticalAlignment','middle');
end

% OPTIONAL: label node numbers
for i = 1:N
    h = text(x(i), y(i), z(i), num2str(nodes(i,1)), ...
             'HorizontalAlignment','center', 'VerticalAlignment','middle');
end
hold off;


%% 2) DRAW OUTPUT: deformed mesh (scaled by plot_scale)
ind_x = 1:3:3*N;
ind_y = 2:3:3*N;
ind_z = 3:3:3*N;

coords_def = coords;
coords_def(:,1) = coords(:,1) + plot_scale * u(ind_x);
coords_def(:,2) = coords(:,2) + plot_scale * u(ind_y);
coords_def(:,3) = coords(:,3) + plot_scale * u(ind_z);

figure(2); clf;
tetramesh(tri, coords_def, 'FaceAlpha', 0.15, 'EdgeColor', 'b');
axis equal;
xlim([min_x-b, max_x+b]);
ylim([min_y-b, max_y+b]);
zlim([min_z-b, max_z+b]);
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
title(sprintf('Deformed mesh (scale = %.2g)', plot_scale));


%% 3) Element centroids (for coloring by element quantity)
centroids = zeros(NN,3);
for e = 1:NN
    node_ids = elements(e,2:5);
    centroids(e,:) = mean(coords(node_ids,:), 1);
end


%% 4) Max principal stress sigma1
figure(3); clf;
scatter3(centroids(:,1), centroids(:,2), centroids(:,3), ...
         40, sigma1, 'filled');
axis equal; grid on;
xlabel('x'); ylabel('y'); zlabel('z');
title('\sigma_1 (max principal stress)');
colorbar;


%% 5) Min principal stress sigma3
figure(4); clf;
scatter3(centroids(:,1), centroids(:,2), centroids(:,3), ...
         40, sigma3, 'filled');
axis equal; grid on;
xlabel('x'); ylabel('y'); zlabel('z');
title('\sigma_3 (min principal stress)');
colorbar;


%% 6) Max principal strain epsi1
figure(5); clf;
scatter3(centroids(:,1), centroids(:,2), centroids(:,3), ...
         40, epsi1, 'filled');
axis equal; grid on;
xlabel('x'); ylabel('y'); zlabel('z');
title('\epsilon_1 (max principal strain)');
colorbar;


%% 7) Min principal strain epsi3
figure(6); clf;
scatter3(centroids(:,1), centroids(:,2), centroids(:,3), ...
         40, epsi3, 'filled');
axis equal; grid on;
xlabel('x'); ylabel('y'); zlabel('z');
title('\epsilon_3 (min principal strain)');
colorbar;