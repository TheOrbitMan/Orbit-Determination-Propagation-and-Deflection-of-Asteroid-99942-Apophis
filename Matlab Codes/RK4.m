%-------------------------------------------------------------------------%
% Runge-Kutta 4 order integrator
% 
% REFERENCE:
% Curtis, H. D.(2013) Orbital Mechanics for Engineering Students, 3. Edition, Elsevier.
% Garcia, A. L. (2000) Numerical Methods for Physics, 2. Edition, Prentice-Hall
%
% Author: Þahin Ulaþ KÖPRÜCÜ
%-------------------------------------------------------------------------%
function [finalstate] = RK4(state,dt,mu,jd0)

pos=[state(1) state(2) state(3)];
vel=[state(4) state(5) state(6)];
acc=feval("kepler",pos,mu)+feval("SRP",pos)+feval("third_bodies",pos,jd0);
k1=[vel acc];

state1=state+(1/2)*k1*dt;jd1=jd0+(1/2)*dt;
pos=[state1(1) state1(2) state1(3)];
vel=[state1(4) state1(5) state1(6)];
acc=feval("kepler",pos,mu)+feval("SRP",pos)+feval("third_bodies",pos,jd1);
k2=[vel acc];

state2=state+(1/2)*k2*dt;jd2=jd0+(1/2)*dt;
pos=[state2(1) state2(2) state2(3)];
vel=[state2(4) state2(5) state2(6)];
acc=feval("kepler",pos,mu)+feval("SRP",pos)+feval("third_bodies",pos,jd2);
k3=[vel acc];

state3=state+k3*dt;jd3=jd0+dt;
pos=[state3(1) state3(2) state3(3)];
vel=[state3(4) state3(5) state3(6)];
acc=feval("kepler",pos,mu)+feval("SRP",pos)+feval("third_bodies",pos,jd3);
k4=[vel acc];

finalstate=state+1/6*(k1+2*k2+2*k3+k4)*dt;

end

