clc;
clear;

T1 = readtable('Takahiro.csv');
T2 = readtable('Josh.csv');

X1 = T1{:,1}; Y1 = T1{:,2}; Z1 = T1{:,3};
Ex1 = T1{:,4}; Ey1 = T1{:,5}; Ez1 = T1{:,6};

X2 = T2{:,1}; Y2 = T2{:,2}; Z2 = T2{:,3};
Ex2 = T2{:,4}; Ey2 = T2{:,5}; Ez2 = T2{:,6};

% Visualize Ex, Ey and Ez
[xq,yq] = meshgrid(-10:0.05:10, -10:0.05:10);
Fieldx1 = griddata(X1,Y1,Ex1,xq,yq,'nearest');
Fieldy1 = griddata(X1,Y1,Ey1,xq,yq,'nearest');
Fieldz1 = griddata(X1,Y1,Ez1,xq,yq,'nearest');

Fieldx2 = griddata(X2,Y2,Ex2,xq,yq,'nearest');
Fieldy2 = griddata(X2,Y2,Ey2,xq,yq,'nearest');
Fieldz2 = griddata(X2,Y2,Ez2,xq,yq,'nearest');

xx1 = [max(max(abs(Fieldx1)))  max(max(abs(Fieldy1)))  max(max(abs(Fieldz1)))]; 
maxfield1 = max(xx1) ; 

xx2 = [max(max(abs(Fieldx2)))  max(max(abs(Fieldy2)))  max(max(abs(Fieldz2)))]; 
maxfield2 = max(xx2) ; 

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
%}

%%%% Magnited of Total Field 1 
%{
Etotal1 = sqrt(abs(Fieldx1).^2 + abs(Fieldy1).^2 + abs(Fieldz1).^2) / sqrt(2);
Etotal_db1 = 20*log10(Etotal1 * 1e6);
figure(7) 
pcolor(xq,yq, Etotal_db1);
shading interp
colorbar
title('|Efield| (dB) Takahiro')

%%%% Magnited of Total Field 2
Etotal2 = sqrt(abs(Fieldx2).^2 + abs(Fieldy2).^2 + abs(Fieldz2).^2) / sqrt(2);
Etotal_db2 = 20*log10(Etotal2 * 1e6);
figure(8) 
pcolor(xq,yq, Etotal_db2);
shading interp
colorbar
title('|Efield| (dB) Josh')

%%%% Magnited of Difference Field
Ediff = abs(Etotal_db2 - Etotal_db1);
figure(9) 
pcolor(xq,yq, Ediff);
shading interp
colorbar
title('|Diff| (dB)')
%}

%%%% Magnited of Total Field 1
Etotal1 = sqrt(abs(Fieldx1).^2 + abs(Fieldy1).^2 + abs(Fieldz1).^2);
Etotal_db1 = 20*log10(Etotal1 * 1e6);
figure(7)
pcolor(xq,yq, Etotal1);
shading interp
colorbar
title('|Efield| (dB) Takahiro')

%%%% Magnited of Total Field 2
Etotal2 = sqrt(abs(Fieldx2).^2 + abs(Fieldy2).^2 + abs(Fieldz2).^2);
Etotal_db2 = 20*log10(Etotal2 * 1e6);
figure(8) 
pcolor(xq,yq, Etotal2);
shading interp
colorbar
title('|Efield| (dB) Josh')

%%%% Magnited of Difference Field
Ediff = abs(Etotal2 - Etotal1);
figure(9) 
pcolor(xq,yq, Ediff);
shading interp
colorbar
title('|Diff| (dB)')
