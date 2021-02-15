%-------------------------------------------------------------------------%
% Computes Sun centered inertial frame to NTW body frame transformation
%
% REFERENCE:
% Vallado, D. A.(2013) Fundamentals of Astrodynamics and Applications, 4. Edition, Microcosm Press.
% 
% Author: Þahin Ulaþ KÖPRÜCÜ
%-------------------------------------------------------------------------%
function [QX] = SCI2NTW(X,Y,Z,VX,VY,VZ)

R=[X Y Z];
V=[VX VY VZ];

T=V/norm(V);
W=cross(R,V)/norm(cross(R,V));
N=cross(T,W);

QX=[N;T;W];
end

