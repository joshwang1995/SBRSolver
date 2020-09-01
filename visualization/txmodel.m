close all; 
clear all ; 

%
% Generate a collection of points to sample the field
% on a constant z-plane (parallel to floor)

xs =  2:.05: 6 ; 
Nx = length(xs) ; 
ys =  2:.05: 6 ; 
Ny = length(ys) ; 

zs = 1 ; 

% Frequency

freq = 2.4e9 ; 

k0 = 2*pi*freq/(3e8) ; % free space wavenumber 


Npts = 0 ; % number of field sampling points

for nx = 1:Nx 
    
    xp = xs(nx)  ; 
    
    for ny = 1:Ny 
        
      Npts = Npts + 1; 
        
      yp = ys(ny) ; 
      
      samplecoord(Npts,1) = xp ; 
      samplecoord(Npts,2) = yp ;
      samplecoord(Npts,3) = zs ; 
      
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
% Transmitter models
%

%
% Step 1: Define transmitter coordinate system
%

% center coordinates
xtc = 4 ; 
ytc = 4 ; 
ztc = 3 ; 

% Unit vectors 
 ux = [0 0 -1]' ;
 uy = [1 0  0]' ; 
 uz = [0 -1 0]' ;


% ux = [1 0 0 ]' ;
% uy = [0 1 0]' ; 
% uz = [0 0 1]' ;

% 
%
% Transformation matrices
%

RTG=[ux(1) uy(1) uz(1); ux(2) uy(2) uz(2); ux(3) uy(3) uz(3)]; 

RGT=RTG' ; 

for n = 1: Npts 
    
    
coordT = RGT*[samplecoord(n,1)-xtc; samplecoord(n,2)-ytc; samplecoord(n,3)-ztc] ; 

RT = sqrt(coordT(1)^2 + coordT(2)^2 + coordT(3)^2) ; 
thetaT = atan(sqrt(coordT(1)^2 + coordT(2)^2)/coordT(3)); 
phiT = atan(coordT(2)/coordT(1)) ; 

st = sin(thetaT); 
sp = sin(phiT) ; 
ct = cos(thetaT); 
cp = cos(phiT) ; 

Rsc = [st*cp ct*cp -sp; st*sp ct*sp cp; ct -st 0] ;

% Field of Hertz dipole 

% ER=0; 
% Ephi = 0 ; 
% Etheta = sin(thetaT)*exp(-j*k0*RT)/RT ; 

% Field of half-wave dipole 

Pt = 1 ; % transmit power (W)

ER=0; 
Ephi = 0; 
Etheta = sqrt(60*Pt)*(cos(pi*cos(thetaT)/2)/sin(thetaT))*exp(-j*k0*RT)/RT ; 


% Field in Tx cartesian coordinates

EcTx = Rsc*[ER; Etheta; Ephi] ; 


% Field in global cartesian coordinates
Ecsample = RTG*EcTx ; 


Ex(n) = Ecsample(1) ; 
Ey(n) = Ecsample(2) ; 
Ez(n) = Ecsample(3) ; 

end

% Visualization

np = 0 ; 

for nx = 1: Nx 
    
    
    for ny = 1: Ny
        
        np = np+1 ; 

      Fieldx (nx,ny) = Ex(np) ; 
      
      Fieldy (nx,ny) = Ey(np) ; 
      
      Fieldz (nx,ny) = Ez(np) ; 
        
      
    end 
    
end

xx = [max(max(abs(Fieldx)))  max(max(abs(Fieldy)))  max(max(abs(Fieldz)))]; 
maxfield = max(xx) ; 

%maxfield = 1 % no normalization 



%%%% NORMALIZED MAGNITUDE PLOTS
figure(1) 

pcolor(abs(Fieldx)/maxfield)
shading interp
caxis([0 1])

title('E_x')

figure(2) 

pcolor(abs(Fieldy)/maxfield)
shading interp
caxis([0 1])
title('E_y')
colorbar

figure(3) 

pcolor(abs(Fieldz)/maxfield)
shading interp
caxis([0 1])
title('E_z')



%%%% PHASE PLOTS
figure(4) 

pcolor(angle(Fieldx)*180/pi)
shading interp
colorbar

title('E_x')

figure(5) 

pcolor(angle(Fieldy)*180/pi)
shading interp
colorbar
title('E_y')

figure(6) 

pcolor(angle(Fieldz)*180/pi)
shading interp
colorbar
title('E_z')








      
      
       




      
      
       
         












