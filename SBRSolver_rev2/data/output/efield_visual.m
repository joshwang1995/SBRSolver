clc;
clear;

T = readtable('ElectricField.csv');

X = T{:,1};
Y = T{:,2};
Z = T{:,3};

Ex = T{:,4};
Ey = T{:,5};
Ez = T{:,6};

Efield = [Ex'; Ey'; Ez'];

% Visualize Ex, Ey and Ez
[xq,yq] = meshgrid(-10:0.05:10, -10:0.05:10);
Fieldx = griddata(X,Y,Ex,xq,yq,'nearest');
Fieldy = griddata(X,Y,Ey,xq,yq,'nearest');
Fieldz = griddata(X,Y,Ez,xq,yq,'nearest');

xx = [max(max(abs(Fieldx)))  max(max(abs(Fieldy)))  max(max(abs(Fieldz)))]; 
maxfield = max(xx) ; 

%maxfield = 1 % no normalization 

%%%% NORMALIZED MAGNITUDE PLOTS
%{
figure(1)
pcolor(xq,yq,abs(Fieldx)/maxfield);
shading interp
caxis([0 1])
title('|E_x|')

figure(2) 
pcolor(xq,yq,abs(Fieldy)/maxfield);
shading interp
caxis([0 1])
title('|E_y|')
colorbar

figure(3) 
pcolor(xq,yq,abs(Fieldz)/maxfield);
shading interp
caxis([0 1])
title('|E_z|')

%%%% PHASE PLOTS
figure(4)
pcolor(xq,yq,angle(Fieldx)*180/pi);
shading interp
colorbar
title('Phase E_x')

figure(5)
pcolor(xq,yq,angle(Fieldy)*180/pi);
shading interp
colorbar
title('Phase E_y')

figure(6) 
pcolor(xq,yq,angle(Fieldz)*180/pi);
shading interp
colorbar
title('Phase E_z')

figure(6) 
pcolor(xq,yq,angle(Fieldz)*180/pi);
shading interp
colorbar
title('Phase E_z')
%}

%%%% Magnited of Total Field
Etotal = sqrt(abs(Fieldx).^2 + abs(Fieldy).^2 + abs(Fieldz).^2) / sqrt(2);
Etotal_db = 20*log10(Etotal * 1e6);
figure(6) 
s = pcolor(xq,yq, Etotal_db);
shading interp
colorbar
title('|E_total|')