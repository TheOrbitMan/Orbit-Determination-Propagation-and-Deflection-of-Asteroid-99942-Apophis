%-------------------------------------------------------------------------%
% Computes acceleration due to Solar Radiation Pressure
% 
% REFERENCE:
% Curtis, H. D.(2013) Orbital Mechanics for Engineering Students, 3. Edition, Elsevier.
%
% Author: Þahin Ulaþ KÖPRÜCÜ
%-------------------------------------------------------------------------%
function [rate] = SRP(pos)

%AU(km)
%au=149597871;
au=149597870.700;

%Reflectivity of the asteroid
cr=1.5;

%Mass of the asteroid(kg)
m=27*(10^9);

%Radius of the asteroid(au)
r=185/(au*1000);

%Absorbing area(au2)
A=pi*r*r;

%Speed of light(au/day)
c=(3*(10^8))*(86400/(au*1000));

%Radius of photosphere(km)
r_photos=696000;

%Sun radiated power intensity at surface(W/au2)
s0=((63.15*(10^6))*(au*1000)^2)*(1/(au*1000))*(86400)*(1/(au*1000))*86400*86400;

%Sun radiated power intensity at asteroid distance(W/au2)
s=s0*(r_photos/(norm(pos)*au))^2;

accel=(((s/c)*cr*A)/m)*(pos/norm(pos));

rate=[accel(1) accel(2) accel(3)];
end

