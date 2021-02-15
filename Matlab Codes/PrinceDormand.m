%-------------------------------------------------------------------------%
% Prince-Dormand integrator
%
% REFERENCE:
% https://core.ac.uk/download/pdf/81989096.pdf
% Garcia, A. L. (2000) Numerical Methods for Physics, 2. Edition, Prentice-Hall
%
% Author: Þahin Ulaþ KÖPRÜCÜ
%-------------------------------------------------------------------------%
function [statefinal] = PrinceDormand(state,dt,mu,jd0)

pos=[state(1) state(2) state(3)];
vel=[state(4) state(5) state(6)];
accel=feval("kepler",pos,mu)+feval("third_bodies",pos,jd0)+feval("SRP",pos);
k1=[vel accel];

state1=state+2/9*k1*dt;jd1=jd0+(2/9)*dt;
pos=[state1(1) state1(2) state1(3)];
vel=[state1(4) state1(5) state1(6)];
accel=feval("kepler",pos,mu)+feval("third_bodies",pos,jd1)+feval("SRP",pos);
k2=[vel accel];

state2=state+1/12*k1*dt+1/4*k2*dt;jd2=jd0+(1/3)*dt;
pos=[state2(1) state2(2) state2(3)];
vel=[state2(4) state2(5) state2(6)];
accel=feval("kepler",pos,mu)+feval("third_bodies",pos,jd2)+feval("SRP",pos);
k3=[vel accel];

state3=state+55/324*k1*dt-25/108*k2*dt+50/81*k3*dt;jd3=jd0+(5/9)*dt;
pos=[state3(1) state3(2) state3(3)];
vel=[state3(4) state3(5) state3(6)];
accel=feval("kepler",pos,mu)+feval("third_bodies",pos,jd3)+feval("SRP",pos);
k4=[vel accel];

state4=state+83/330*k1*dt-13/22*k2*dt+61/66*k3*dt+9/110*k4*dt;jd4=jd0+(2/3)*dt;
pos=[state4(1) state4(2) state4(3)];
vel=[state4(4) state4(5) state4(6)];
accel=feval("kepler",pos,mu)+feval("third_bodies",pos,jd4)+feval("SRP",pos);
k5=[vel accel];

state5=state-19/28*k1*dt+9/4*k2*dt+1/7*k3*dt-27/7*k4*dt+22/7*k5*dt;jd5=jd0+dt;
pos=[state5(1) state5(2) state5(3)];
vel=[state5(4) state5(5) state5(6)];
accel=feval("kepler",pos,mu)+feval("third_bodies",pos,jd5)+feval("SRP",pos);
k6=[vel accel];

statefinal=state+(19/200*k1+3/5*k3-243/400*k4+33/40*k5+7/80*k6)*dt;
end

