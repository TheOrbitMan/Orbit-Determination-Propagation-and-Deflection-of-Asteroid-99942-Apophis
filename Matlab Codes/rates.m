%-------------------------------------------------------------------------%
% Rates function
%
% Author: Þahin Ulaþ KÖPRÜCÜ 
%-------------------------------------------------------------------------%
function dfdt = rates(t,f,mu)

r=f(1:3)';
v=f(4:6)';

accel1=kepler(r,mu);
accel2=SRP(r);
accel3=third_bodies(r,t);

dx=v';
d2x=accel1'+accel2'+accel3';

dfdt=[dx;d2x];

end

