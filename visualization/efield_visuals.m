clc;
clear;

T = readtable('efield.csv');
%T = sortrows(T);

X = T{:,1};
Y = T{:,2};
Z = T{:,3};

Ex = complex(T{:,4},T{:,5});
Ey = complex(T{:,6},T{:,7});
Ez = complex(T{:,8},T{:,9});

Efield = [Ex'; Ey'; Ez'];

%{
% Draw the points of intersection on the defined plane 
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
xlabel('x'); ylabel('y'); zlabel('z'); 
%}

% Visualize Ex, Ey and Ez
[xq,yq] = meshgrid(2:.05:6, 2:.05:6);
Fieldx = griddata(X,Y,Ex,xq,yq);
Fieldy = griddata(X,Y,Ey,xq,yq);
Fieldz = griddata(X,Y,Ez,xq,yq);

xx = [max(max(abs(Fieldx)))  max(max(abs(Fieldy)))  max(max(abs(Fieldz)))]; 
maxfield = max(xx) ; 

%maxfield = 1 % no normalization 

%%%% NORMALIZED MAGNITUDE PLOTS
figure(1) 
pcolor(xq,yq,abs(Fieldx)/maxfield);
shading interp
caxis([0 1])
title('E_x')

figure(2) 
pcolor(xq,yq,abs(Fieldy)/maxfield);
shading interp
caxis([0 1])
title('E_y')
colorbar

figure(3) 
pcolor(xq,yq,abs(Fieldz)/maxfield);
shading interp
caxis([0 1])
title('E_z')

%%%% PHASE PLOTS
figure(4)
pcolor(xq,yq,angle(Fieldx)*180/pi);
shading interp
colorbar
title('E_x')

figure(5)
pcolor(xq,yq,angle(Fieldy)*180/pi);
shading interp
colorbar
title('E_y')

figure(6) 
pcolor(xq,yq,angle(Fieldz)*180/pi);
shading interp
colorbar
title('E_z')