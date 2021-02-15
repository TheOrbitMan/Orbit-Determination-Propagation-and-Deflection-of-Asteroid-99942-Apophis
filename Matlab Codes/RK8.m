%-------------------------------------------------------------------------%
% Runge-Kutta 8 order integrator
% 
% REFERENCE:
% https://ntrs.nasa.gov/citations/19760017203
% Garcia, A. L. (2000) Numerical Methods for Physics, 2. Edition, Prentice-Hall
%
% Author: Þahin Ulaþ KÖPRÜCÜ
%-------------------------------------------------------------------------%
function [finalstate] = RK8(state,dt,mu,jd0)

pos=[state(1) state(2) state(3)];
vel=[state(4) state(5) state(6)];
acc=feval("kepler",pos,mu)+feval("SRP",pos)+feval("third_bodies",pos,jd0);
k1=[vel acc];

state1=state+(4/27)*k1*dt;jd1=jd0+(4/27)*dt;
pos=[state1(1) state1(2) state1(3)];
vel=[state1(4) state1(5) state1(6)];
acc=feval("kepler",pos,mu)+feval("SRP",pos)+feval("third_bodies",pos,jd1);
k2=[vel acc];

state2=state+(1/18)*(k1+3*k2)*dt;jd2=jd0+(2/9)*dt;
pos=[state2(1) state2(2) state2(3)];
vel=[state2(4) state2(5) state2(6)];
acc=feval("kepler",pos,mu)+feval("SRP",pos)+feval("third_bodies",pos,jd2);
k3=[vel acc];

state3=state+(1/12)*(k1+3*k3)*dt;jd3=jd0+(1/3)*dt;
pos=[state3(1) state3(2) state3(3)];
vel=[state3(4) state3(5) state3(6)];
acc=feval("kepler",pos,mu)+feval("SRP",pos)+feval("third_bodies",pos,jd3);
k4=[vel acc];

state4=state+(1/8)*(k1+3*k4)*dt;jd4=jd0+(1/2)*dt;
pos=[state4(1) state4(2) state4(3)];
vel=[state4(4) state4(5) state4(6)];
acc=feval("kepler",pos,mu)+feval("SRP",pos)+feval("third_bodies",pos,jd4);
k5=[vel acc];

state5=state+(1/54)*(13*k1-27*k3+42*k4+8*k5)*dt;jd5=jd0+(2/5)*dt;
pos=[state5(1) state5(2) state5(3)];
vel=[state5(4) state5(5) state5(6)];
acc=feval("kepler",pos,mu)+feval("SRP",pos)+feval("third_bodies",pos,jd5);
k6=[vel acc];

state6=state+(1/4320)*(389*k1-54*k3+966*k4-824*k5+243*k6)*dt;jd6=jd0+(1/6)*dt;
pos=[state6(1) state6(2) state6(3)];
vel=[state6(4) state6(5) state6(6)];
acc=feval("kepler",pos,mu)+feval("SRP",pos)+feval("third_bodies",pos,jd6);
k7=[vel acc];

state7=state+(1/20)*(-234*k1+81*k3-1164*k4+656*k5-122*k6+800*k7)*dt;jd7=jd0+(1)*dt;
pos=[state7(1) state7(2) state7(3)];
vel=[state7(4) state7(5) state7(6)];
acc=feval("kepler",pos,mu)+feval("SRP",pos)+feval("third_bodies",pos,jd7);
k8=[vel acc];

state8=state+(1/288)*(-127*k1+18*k3-678*k4+456*k5-9*k6+576*k7+4*k8)*dt;jd8=jd0+(5/6)*dt;
pos=[state8(1) state8(2) state8(3)];
vel=[state8(4) state8(5) state8(6)];
acc=feval("kepler",pos,mu)+feval("SRP",pos)+feval("third_bodies",pos,jd8);
k9=[vel acc];

state9=state+(1/820)*(1481*k1-81*k3+7104*k4-3376*k5+72*k6-5040*k7-60*k8+720*k9)*dt;jd9=jd0+(1)*dt;
pos=[state9(1) state9(2) state9(3)];
vel=[state9(4) state9(5) state9(6)];
acc=feval("kepler",pos,mu)+feval("SRP",pos)+feval("third_bodies",pos,jd9);
k10=[vel acc];

finalstate=state+(1/840)*(41*k1+27*k4+272*k5+27*k6+216*k7+216*k9+41*k10)*dt;

end

