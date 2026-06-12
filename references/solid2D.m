function [u,epsi,sigma] = solid2D(nodes,elements,load1,load2,bcs)
%NUMBER of NODES
N = length(nodes(:,1));

%X- og y- cordinates of noeds
x = nodes(:,2);
y = nodes(:,3);

%Number of elements
n = length(elements(:,1));

%Global stiffness matrix
K = zeros(2*N,2*N);
k1 = 0;
save('stiffness_matrices','k1')
for i = 1:n,
    E = elements(i,5);
    nyy = elements(i,6);
    h = elements(i,7);
    
    %Nodes
    I = elements(i,2);
    J = elements(i,3);
    Q = elements(i,4);
    
    x1 = nodes(I,2);
    x2 = nodes(J,2);
    x3 = nodes(Q,2);
    y1 = nodes(I,3);
    y2 = nodes(J,3);
    y3 = nodes(Q,3);
    b1 = y2-y3; b2 = y3-y1; b3 = y1-y2;
    c1 = x3-x2; c2 = x1-x3; c3 = x2-x1;
    f1 = x2*y3-x3*y2; f2 = x3*y1-x1*y3; f3 = x1*y2-x2*y1;
    A = (f1+f2+f3)/2;
    
    B = 1/(2*A)*[b1 0 b2 0 b3 0; 0 c1 0 c2 0 c3; c1 b1 c2 b2 c3 b3]';
    C = E/(1-nyy^2)*[1 nyy 0; nyy 1 0; 0 0 (1-nyy)/2]
    k = h*A*B*C*B'
    
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
    I1 = I*2-1;
    I2 = I*2;
    J1 = J*2-1;
    J2 = J*2;
    Q1 = Q*2-1;
    Q2 = Q*2;
    index = [I1 I2 J1 J2 Q1 Q2];
    
    %Global stiffness matrix
    K_fake = zeros(2*N,2*N);
    K_fake(index,index) = k;
    K = K+K_fake;
end


%Load vector
NN = length(load1(:,1));
rp = zeros(2*N,1);
for i = 1:NN,
    f = load1(i,2:3);
    I = 2*load1(i,1) - 1;
    J = 2*load1(i,1);
    index = [I J];
    rp(index) = f;
end

%Load - distributed 
%The load is specified in terms of element edge numbers. They are numbered
%1, 2 and 3

if isempty(load2) == 0
NN = length(load2(:,1));
rq = zeros(2*N,1);
q_teikna = zeros(NN,8);
for i = 1:NN,
    %Determine which element is loaded
    nr_e = load2(i,1);
    
    %Determine which edge is loaded
    nr_side = load2(i,4);
    
    %Determine the node numbers that are affected
    ii = elements(nr_e,2);
    jj = elements(nr_e,3);
    qq = elements(nr_e,4);
 
    
    %Determine which DOFs are affected
    I1 = 2*ii - 1;
    I2 = 2*ii;
    J1 = 2*jj - 1;
    J2 = 2*jj;
    Q1 = 2*qq-1;
    Q2 = 2*qq;

    %Load is qn=normal stress or qt=shear stress
    qn = load2(i,2);
    qt = load2(i,3);
    
    %nodes
    x1 = nodes(ii,2);
    x2 = nodes(jj,2);
    x3 = nodes(qq,2);
    y1 = nodes(ii,3);
    y2 = nodes(jj,3);
    y3 = nodes(qq,3);
    
    %Edge lengths
    L12 = sqrt((x1-x2)^2+(y1-y2)^2);
    L23 = sqrt((x2-x3)^2+(y2-y3)^2);
    L31 = sqrt((x3-x1)^2+(y3-y1)^2);
    
    if nr_side == 1
        nx = (y2-y1)/L12; ny = -(x2-x1)/L12;
        rq1 = h*L12/2*[nx*qn-ny*qt ny*qn+nx*qt nx*qn-ny*qt ny*qn+nx*qt 0 0];
        input = [rq1]';
        q_teikna(i,:) = [x1 x2 y1 y2 nx ny qn qt];
    end
    
    if nr_side == 2
        nx = (y3-y2)/L23;
        ny = -(x3-x2)/L23;
        rq2 = h*L23/2*[0 0 nx*qn-ny*qt ny*qn+nx*qt nx*qn-ny*qt ny*qn+nx*qt];
        input = [rq2]';
        q_teikna(i,:) = [x2 x3 y2 y3 nx ny qn qt];
    end
    
    if nr_side == 3
        nx = (y1-y3)/L31; ny = -(x1-x3)/L31;
        rq3 = h*L31/2*[nx*qn-ny*qt ny*qn+nx*qt 0 0 nx*qn-ny*qt ny*qn+nx*qt];
        input = [rq3]';
        q_teikna(i,:) = [x1 x3 y1 y3 nx ny qn qt];
    end
    
    index = [I1 I2 J1 J2 Q1 Q2];
    rq(index) = rq(index) + input; 
end

%Adding up load vectors
F = rp + rq;
else
F = rp;    
end

%BCs
NN = length(bcs(:,1));
F_old = F;
for i = 1:NN,
    r = bcs(i,2:3);
    I = 2*bcs(i,1) - 1;
    J = 2*bcs(i,1);
    index = [I J];
    F(index) = r;

    %The effect of BCs on the global stiffness matrix
    if isnan(r(1)) == 0,
        K(I,:) = zeros(1,2*N);
        K(I,I) = 1;
    end

    if isnan(r(2)) == 0,
        K(J,:) = zeros(1,2*N);
        K(J,J) = 1;
    end
end
I = find(isnan(F));
F(I) = F_old(I);
F = F
%pause

%Solution
u = K\F;


%Stress and strain
epsi = zeros(n,4);
epsi(:,1) = 1:n;
sigma = zeros(n,4);
sigma(:,1) = 1:n;

for i = 1:n,
    E = elements(i,5);
    nyy = elements(i,6);
    h = elements(i,7);
    
    %Nodes
    I = elements(i,2);
    J = elements(i,3);
    Q = elements(i,4);
    
    I1 = 2*I - 1;
    I2 = 2*I;
    J1 = 2*J - 1;
    J2 = 2*J;
    Q1 = 2*Q-1;
    Q2 = 2*Q;
    index = [I1 I2 J1 J2 Q1 Q2];
    
    x1 = nodes(I,2);
    x2 = nodes(J,2);
    x3 = nodes(Q,2);
    y1 = nodes(I,3);
    y2 = nodes(J,3);
    y3 = nodes(Q,3);
    b1 = y2-y3; b2 = y3-y1; b3 = y1-y2;
    c1 = x3-x2; c2 = x1-x3; c3 = x2-x1;
    f1 = x2*y3-x3*y2; f2 = x3*y1-x1*y3; f3 = x1*y2-x2*y1;
    A = (f1+f2+f3)/2;
    
    B = 1/(2*A)*[b1 0 b2 0 b3 0; 0 c1 0 c2 0 c3; c1 b1 c2 b2 c3 b3]';
    C = E/(1-nyy^2)*[1 nyy 0; nyy 1 0; 0 0 (1-nyy)/2];
    k = h*A*B*C*B';
    
    d = u(index);
    epsi(i,2:4) = B'*d;
    sigma(i,2:4) = C*epsi(i,2:4)';
end


%DRAW DRAW
n = length(elements(:,1));
min_x = min(min(x));
max_x = max(max(x));
min_y = min(min(y));
max_y = max(max(y));
b = 0.1*max(mean([min_x max_x]),mean([min_y max_y]));

figure(1)
clf
set(gca,'xlim',[min_x-b max_x+b])

%DRAW input
figure(1)
hold on
tri = elements(:,2:4);
h = triplot(tri,x,y,'k');
%set(h,'markeredgecolor',[0 0 0])
xlabel('x')
ylabel('y')
NN = length(elements(:,1));
for i =1:NN,
    I = elements(i,2);
    J = elements(i,3);
    Q = elements(i,4);
    
    x1 = nodes(I,2);
    x2 = nodes(J,2);
    x3 = nodes(Q,2);
    y1 = nodes(I,3);
    y2 = nodes(J,3);
    y3 = nodes(Q,3);
    ave_x = mean([x1 x2 x3]);
    ave_y = mean([y1 y2 y3]);
    h = text(ave_x,ave_y,num2str(i));
    set(h,'HorizontalAlignment','center')
    set(h,'VerticalAlignment','middle')
end

for i = 1:N,
    xx = nodes(i,2);
    yy = nodes(i,3);
    h = text(xx,yy,num2str(nodes(i,1)));
    set(h,'HorizontalAlignment','center')
    set(h,'VerticalAlignment','middle')
end

hold off


%DRAW BCs
figure(1)
hold on
NN = length(bcs(:,1));
for i = 1:NN,
    nrI = bcs(i,1);
    r = bcs(i,2:3);
    if isnan(r(1)) == 0
        x1 = nodes(nrI,2);
        y1 = nodes(nrI,3);
        d = b/2;
        h = plot([x1 x1-2*d/3 x1+2*d/3 x1],[y1 y1-d y1-d y1],'r');
        set(h,'linewidth',[1.25])
    end
    if isnan(r(2)) == 0
        x1 = nodes(nrI,2);
        y1 = nodes(nrI,3);
        d = b/2;
        h = plot([x1 x1-d x1-d x1],[y1 y1-2*d/3 y1+2*d/3 y1],'r');
        set(h,'linewidth',[1.25])
    end
end
hold off


%DRAW force input
figure(1)
hold on
NN = length(load1(:,1));
for i = 1:NN,
    nrI = load1(i,1);
    r = load1(i,2:3);
    
    %Horizontal force
    if r(1) ~= 0
        x1 = nodes(nrI,2);
        y1 = nodes(nrI,3);
        if r(1) < 0
            h = plot([x1 x1-b],[y1 y1],'r');
            set(h,'linewidth',[1.5])
            h = plot([x1-b],[y1],'r<');
            set(h,'markerfacecolor',[1 0 0])
            set(h,'markersize',[8])
        end
        if r(1) > 0
            h = plot([x1 x1+b],[y1 y1],'r');
            set(h,'linewidth',[1.5])
            h = plot([x1+b],[y1],'r>');
            set(h,'markerfacecolor',[1 0 0])
            set(h,'markersize',[8])
        end
    end
    
    %Vertical force
    if r(2) ~= 0
        if r(2) < 0
            x1 = nodes(nrI,2);
            h = plot([x1 x1],[y1 y1-b],'r');
            set(h,'linewidth',[1.5])
            h = plot([x1],[y1-b],'rv');
            set(h,'markerfacecolor',[1 0 0])
            set(h,'markersize',[8])
        end
        if r(2) > 0
            x1 = nodes(nrI,2);
            h = plot([x1 x1],[y1 y1+b],'r');
            set(h,'linewidth',[1.5])
            h = plot([x1],[y1+b],'r^');
            set(h,'markerfacecolor',[1 0 0])
            set(h,'markersize',[8])
        end
    end
end
hold off

%Draw distributed load
figure(1)
hold on
axis equal
if isempty(load2) == 0
NNN = length(load2(:,1));
for i = 1:NNN,
    x1 = q_teikna(i,1);
    x2 = q_teikna(i,2);
    y1 = q_teikna(i,3);
    y2 = q_teikna(i,4);
    nx = q_teikna(i,5);
    ny = q_teikna(i,6);
    qn = q_teikna(i,7);
    qt = q_teikna(i,8)
    QQ = 10;
    if qn ~= 0
        h = quiver(linspace(x1,x2,QQ),linspace(y1,y2,QQ),sign(qn)*nx*linspace(1,1,QQ),sign(qn)*ny*linspace(1,1,QQ),0.15);
        set(h,'color',[1 0 1])
        set(h,'linewidth',[1.5])
        set(h,'markersize',[10])
        set(h,'markerfacecolor',[1 0 1])
    end
    if qt ~= 0
        h = quiver(linspace(x1,x2,QQ),linspace(y1,y2,QQ),-sign(qt)*ny*linspace(1,1,QQ),sign(qt)*nx*linspace(1,1,QQ),0.15);
        set(h,'color',[0 0 0])
        set(h,'linewidth',[1.5])
        set(h,'markersize',[10])
        set(h,'markerfacecolor',[0 0 0])
    end
end
end
print -dpng figure1
hold off

%DRAW normal stress in x-direction
figure(2)
clf
axis equal
hold on
tri = elements(:,2:4);
h = trisurf(tri,x,y,ones(size(x)),sigma(:,2));
%for i = 1:length(sigma(:,1))
%    avex = mean(x(tri(i,:)));
%    avey = mean(y(tri(i,:)));
%    h = text(avex,avey,num2str(sigma(i,2)));
%    set(h,'Color',[1 1 1])
%    set(h,'fontsize',14)
%end
colorbar('vert')
colormap jet
%set(h,'markeredgecolor',[0 0 0])
h = title('\sigma_x');
set(h,'fontsize',16)
xlabel('x')
ylabel('y')
max_sigma = max(max(abs(sigma(:,2))))
print -dpng figure2
hold off

%DRAW normal stress in y-direction
figure(3)
clf
axis equal
hold on
tri = elements(:,2:4);
h = trisurf(tri,x,y,ones(size(x)),sigma(:,3));
colorbar('vert')
colormap jet
%set(h,'markeredgecolor',[0 0 0])
h = title('\sigma_y');
set(h,'fontsize',16)
xlabel('x')
ylabel('y')
print -dpng figure3
hold off


%DRAW shear stress in xy-plane
figure(4)
clf
axis equal
hold on
tri = elements(:,2:4);
h = trisurf(tri,x,y,ones(size(x)),sigma(:,4));
colorbar('vert')
colormap jet
%set(h,'markeredgecolor',[0 0 0])
h = title('\tau_{xy}');
set(h,'fontsize',16)
xlabel('x')
ylabel('y')
print -dpng figure4
hold off

%DRAW normal strain in x-direction
figure(5)
clf
axis equal
hold on
tri = elements(:,2:4);
h = trisurf(tri,x,y,ones(size(x)),epsi(:,2));
colorbar('vert')
colormap jet
%set(h,'markeredgecolor',[0 0 0])
h = title('\epsilon_x');
set(h,'fontsize',16)
xlabel('x')
ylabel('y')
print -dpng figure5
hold off

%DRAW normal strain in y-direction
figure(6)
clf
axis equal
hold on
tri = elements(:,2:4);
h = trisurf(tri,x,y,ones(size(x)),epsi(:,3));
colorbar('vert')
colormap jet
%set(h,'markeredgecolor',[0 0 0])
h = title('\epsilon_y');
set(h,'fontsize',16)
xlabel('x')
ylabel('y')
print -dpng figure6
hold off

%DRAW shear strain in xy-plane
figure(7)
clf
axis equal
hold on
tri = elements(:,2:4);
h = trisurf(tri,x,y,ones(size(x)),epsi(:,4));
colorbar('vert')
colormap jet
%set(h,'markeredgecolor',[0 0 0])
h = title('\gamma_{xy}');
set(h,'fontsize',16)
xlabel('x')
ylabel('y')
print -dpng figure7
hold off





if 1 == 10000000
for i = 1:n,
    qn = load2(i,2);
    qt = load2(i,3);
    
    %Finn fyrst frelsigráður elementsinnar
    I = elements(i,2);
    I1 = I*2-1;
    I2 = I*2;
    J = elements(i,3);
    J1 = J*2-1;
    J2 = J*2;
    
    %nodes elementsinnar
    x1 = nodes(I,2);
    x2 = nodes(J,2);
    y1 = nodes(I,3);
    y2 = nodes(J,3);
    L = sqrt((x1-x2)^2 + (y1-y2)^2);
    
    %Pilla út frelsigráðulausnir
    index = [I1 I2 J1 J2];
    d = u(index);
    
    %Teikna jafndreift álag
    figure(1)
    hold on
    QQ = 10;
    quiver(linspace(x1,x2,QQ),linspace(y1,y2,QQ),linspace(qn,qn,QQ),linspace(qt,qt,QQ))
    hold off
    
    
    if 1 == 1000000
    %Reikna út vægi
    figure(2)
    hold on
    h = plot([x1 x2],[0 0],'k');
    set(h,'linewidth',[1.25])
    h = plot([x1 x2],[0 0],'ko');
    set(h,'markersize',[6])
    set(h,'markerfacecolor','k')
    h = text(x1,0,num2str(I),'verticalalignment','bottom','horizontalalignment','center');
    h = text(x2,0,num2str(J),'verticalalignment','bottom','horizontalalignment','center');
   
    %text(x,y,sprintf('%4.1f',f(i)),'Horizontalalignment','center')
    E = elements(i,4);
    Iz = elements(i,5);
    x = linspace(x1,x2,100);
    s = x-x1;
    
    %Jafndreift álag
    M = E*Iz*((-6/L^2+12*s/L^3)*d(1) + (-4/L+6*s/L^2)*d(2) ...
        + (6/L^2-12*s/L^3)*d(3) + (-2/L+6*s/L^2)*d(4)) ...
        + q_vec(i)/24*(2*L^2 - 12*L*s + 12*s.^2);
    
    h = plot(x,M/1000,'r');
    set(h,'linewidth',[1.25])
    ylabel('M [kNm]')
    xlabel('x [m]')
    grid on
    hold off
    
    %Reikna út niðurbeygju
    figure(3)
    hold on
    h = plot([x1 x2],[0 0],'k');
    set(h,'linewidth',[1.25])
    h = plot([x1 x2],[0 0],'ko');
    set(h,'markersize',[6])
    set(h,'markerfacecolor','k')
    h = text(x1,0,num2str(I),'verticalalignment','bottom','horizontalalignment','center');
    h = text(x2,0,num2str(J),'verticalalignment','bottom','horizontalalignment','center');

    v = (2*s.^3/L^3-3*s.^2/L^2+1)*d(1) + (s.^3/L^2-2*s.^2/L+s)*d(2) ...
        + (3*s.^2/L^2-2*s.^3/L^3)*d(3) + (s.^3/L^2-s.^2/L)*d(4) ...
        + q_vec(i)/(24*E*Iz)*(L^2*s.^2 - 2*L*s.^3 + s.^4);
    h = plot(x,v*1000,'r');
    set(h,'linewidth',[1.25])
    ylabel('v [mm]')
    xlabel('x [m]')
    grid on
    hold off
    end
end





ind_odd = 1:2:2*N;
    ind_slett = 2:2:2*N;
    X = nodes(:,2);
    Y = nodes(:,3);
    d = 0.15;
    Fx = F(ind_odd);
    Fx = Fx./abs(Fx);
    Fy = F(ind_slett);
    Fy = Fy./abs(Fy);

    hold on
    %quiver(X,Y,Fx,zeros(N,1),d)
    %quiver(X,Y,zeros(N,1),Fy,d)
    quiver(X,Y,Fx,zeros(N,1),d)
    quiver(X,Y,zeros(N,1),Fy,d)
    hold off





%Teikna formbreytt virki
figure(2)
clf
[tt] = teikna(nodes,elements,[1 0 0],'',0,0,zeros(n,1));

hold on
delta_nodes = zeros(N,2);
delta_nodes(:,1) = u(ind_odd);
delta_nodes(:,2) = u(ind_slett);
nodes(:,2:3) = nodes(:,2:3) + skolun*delta_nodes;
[tt] = teikna(nodes,elements,[0 0 1],'Formbreytt virki',1,1,zeros(n,1));
hold off

%Teikna stangarkrafta
figure(3)
clf
hold on
[tt] = teikna(nodes,elements,[0 0 1],'Stangarkraftar [kN]',1,0,f/1000);
hold off

L = L'
end