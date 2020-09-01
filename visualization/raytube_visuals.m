clc;
clear;

T = readtable('raytube.csv');

x0 = T{:,1};
y0 = T{:,2};
z0 = T{:,3};
xr = T{:,4};
yr = T{:,5};
zr = T{:,6};
ray_tube = T{:,7};

t = 1;
xp = x0+t*xr;
yp = y0+t*yr;
zp = z0+t*zr;

[unique_groups, ~, group_idx] = unique(ray_tube);
num_groups = size(unique_groups, 1);
group_names = cellstr( num2str( unique_groups ) );
cmap = jet(num_groups);    %or build a custom color map

hold on;
for i=1:num_groups
    first = find(ray_tube == i, 1, 'first');
    last = find(ray_tube == i, 1, 'last');
    if i<=12
        quiver3(x0(first:last,1),y0(first:last,1),z0(first:last,1),xr(first:last,1),yr(first:last,1),zr(first:last,1),'color',[0 0 0]);
    else
        quiver3(x0(first:last,1),y0(first:last,1),z0(first:last,1),xr(first:last,1),yr(first:last,1),zr(first:last,1),'color',cmap(i,:));
    end
    
end
%pointsize = 50;
%scatter3(xr, yr, zr, pointsize, group_idx);

%colormap( cmap );
%quiver3(x0,y0,z0,xr,yr,zr,'color',ray_tube);
%xlabel('x'); ylabel('y'); zlabel('z'); 
