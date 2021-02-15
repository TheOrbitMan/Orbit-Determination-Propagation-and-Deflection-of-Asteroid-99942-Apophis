%-------------------------------------------------------------------------%
% Computes Ecliptic frame to Equatorial frame transformation
%
% REFERENCE:
% Boulet, D. L.(1991) Methods of Orbit Determination for the Microcomputer, 1. Edition, Willmann-Bell.
% 
% Author: Þahin Ulaþ KÖPRÜCÜ
%-------------------------------------------------------------------------%
function [req,veq] = Ec2Eq(r,v,jd)

T=(jd - 2451545)/36525;
E=23.439291-0.0130042*T-0.00000016*T*T;

x=r(1);
y=r(2)*cosd(E)-r(3)*sind(E);
z=r(2)*sind(E)+r(3)*cosd(E);

xdot=v(1);
ydot=v(2)*cosd(E)-v(3)*sind(E);
zdot=v(2)*sind(E)+v(3)*cosd(E);

veq=[xdot ydot zdot];
req=[x y z];

end

