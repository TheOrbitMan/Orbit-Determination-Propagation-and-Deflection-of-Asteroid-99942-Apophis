%-------------------------------------------------------------------------%
% Laplace's orbit determination
%
% REFERENCE:
% Vallado, D. A.(2013) Fundamentals of Astrodynamics and Applications, 4. Edition, Microcosm Press.
%
% Author: Þahin Ulaþ KÖPRÜCÜ
%-------------------------------------------------------------------------%
function [r2,v2] = LAPLACE(rhohat1, rhohat2, rhohat3, R1, R2, R3, jd1, jd2, jd3, mu,xi)

k=1;
jd1=jd1*k;jd2=jd2*k;jd3=jd3*k;

rhohat2dot=(2*jd2-jd2-jd3)/((jd1-jd2)*(jd1-jd3))*rhohat1+...
    (2*jd2-jd1-jd3)/((jd2-jd1)*(jd2-jd3))*rhohat2+...
    (2*jd2-jd1-jd2)/((jd3-jd1)*(jd3-jd2))*rhohat3;

rhohat2ddot=(2)/((jd1-jd2)*(jd1-jd3))*rhohat1+...
    (2)/((jd2-jd1)*(jd2-jd3))*rhohat2+...
    (2)/((jd3-jd1)*(jd3-jd2))*rhohat3;

R2dot=(2*jd2-jd2-jd3)/((jd1-jd2)*(jd1-jd3))*R1+...
    (2*jd2-jd1-jd3)/((jd2-jd1)*(jd2-jd3))*R2+...
    (2*jd2-jd1-jd2)/((jd3-jd1)*(jd3-jd2))*R3;

R2ddot=(2)/((jd1-jd2)*(jd1-jd3))*R1+...
    (2)/((jd2-jd1)*(jd2-jd3))*R2+...
    (2)/((jd3-jd1)*(jd3-jd2))*R3;

D=det([rhohat2' rhohat2dot' rhohat2ddot'])*2;
D1=det([rhohat2' rhohat2dot' R2ddot']);
D2=det([rhohat2' rhohat2dot' R2']);
C=dot(rhohat2,R2);

a=(4*C*D1)/D-4*(D1^2/D^2)-norm(R2)^2;
b=((4*C*D2)/D-(8*D1*D2)/D^2)*mu;
c=-(4*mu*mu*D2^2)/D^2;

%NEWTON-RAPHSON
error=100;
limit=10e-8;
x=xi;
while error>limit
    x0=x;
    x=x0-(x0^8+a*x0^6+b*x0^3+c)/(8*x0^7+6*a*x0^5+3*b*x0^2);
    error=(abs(x-x0)/x0)*100;
end

r=x;

rho2=-2*(D1/D)-2*(mu/r^3)*(D2/D);

D3=det([rhohat2' R2ddot' rhohat2ddot']);
D4=det([rhohat2' R2' rhohat2ddot']);

rhodot=-D3/D-(mu/r^3)*(D4/D);

r2=rhohat2*rho2+R2;
v2=rhodot*rhohat2+rho2*rhohat2dot+R2dot;
end

