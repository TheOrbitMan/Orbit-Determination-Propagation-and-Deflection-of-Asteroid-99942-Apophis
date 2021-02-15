%-------------------------------------------------------------------------%
% Two-Body Analytical Propagator
%
% Author: Þahin Ulaþ KÖPRÜCÜ 
%-------------------------------------------------------------------------%
function [r,v] = Two_Body_Propagator(r0,v0,t,mu)

coe=coe_from_sv(r0,v0,mu);
TA0=coe(5);                              %rad
e=coe(1);
a=coe(6);                                %au
RA=coe(2);                               %rad
i=coe(3);                                %rad
w=coe(4);                                %rad
T=((2*pi)/sqrt(mu))*a^(1.5);             %day
n=(2*pi)/T;                              %rad/day
h=sqrt(a*mu*(1-e*e));                    %au^2/day

%Eccentric anomaly at t=0 (rad)
E0=2*atan(sqrt((1-e)/(1+e))*tan(TA0/2)); 

%Mean anomaly at t=0 (rad)
M0=E0-e*sin(E0);

%Mean anomaly at t (rad)
M=mod(M0+n*t,2*pi);

%Solve the kepler equation to find eccentric anomaly at t (rad)
E=kepler_E(e, M);

%True anomaly at t (rad)
TA=2*atan(sqrt((1+e)/(1-e))*tan(E/2));

coe=[h e RA i w TA];

[r, v] = sv_from_coe(coe,mu);

end