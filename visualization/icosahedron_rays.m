clc;
clear;

T = readtable('rays.csv');
T2 = readtable('output.csv');

X = T2{:,1};
Y = T2{:,2};
Z = T2{:,3};

x0 = T{:,1};
y0 = T{:,2};
z0 = T{:,3};
xr = T{:,4};
yr = T{:,5};
zr = T{:,6};

t=100;
xp = x0+t*xr;
yp = y0+t*yr;
zp = z0+t*zr;

p1 = [6,2,1];
p2 = [6,6,1];
p3 = [2,6,1];
p4 = [2,2,1];

x_plane = [p1(1) p2(1) p3(1) p4(1)];
y_plane = [p1(2) p2(2) p3(2) p4(2)];
z_plane = [p1(3) p2(3) p3(3) p4(3)];

hold on;
fill3(x_plane, y_plane, z_plane, 1);
scatter3(X,Y,Z);
quiver3(x0,y0,z0,xr,yr,zr);
xlabel('x'); ylabel('y'); zlabel('z'); 

%scatter3(xp,yp,zp);

