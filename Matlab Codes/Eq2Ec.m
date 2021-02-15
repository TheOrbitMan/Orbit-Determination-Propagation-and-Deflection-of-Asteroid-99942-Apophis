%-------------------------------------------------------------------------%
% Computes Equatorial frame to Ecliptic frame transformation
%
% REFERENCE:
% Boulet, D. L.(1991) Methods of Orbit Determination for the Microcomputer, 1. Edition, Willmann-Bell.
% 
% Author: Þahin Ulaþ KÖPRÜCÜ
%-------------------------------------------------------------------------%
function [rec,vec] = Eq2Ec(r,v,jd)

T=(jd - 2451545)/36525;
E=23.439291-0.0130042*T-0.00000016*T*T;

x=r(1);
y=r(3)*sind(E)+r(2)*cosd(E);
z=r(3)*cosd(E)-r(2)*sind(E);

xdot=v(1);
ydot=v(3)*sind(E)+v(2)*cosd(E);
zdot=v(3)*cosd(E)-v(2)*sind(E);

rec=[x y z];
vec=[xdot ydot zdot];
end

