%-------------------------------------------------------------------------%
% Plots the deflection graph from the saved variables
%
% Author: Þahin Ulaþ KÖPRÜCÜ
%-------------------------------------------------------------------------%
clc;clear;
%% delv=0.5 m/s
%
load N_direction05.mat
distance1=distance_TCA;
load T_direction05.mat
distance2=distance_TCA;
load W_direction05.mat
distance3=distance_TCA;
load times.mat

plot(time,distance1)
hold on
plot(time,distance2)
hold on
plot(time,distance3)
xlabel('Deflection times (day)')
ylabel('Distance from Earth at TCA (Re)')
legend('N','T','W')
title('Deflection with 0.5 m/s')
figure
%}
%% delv=1 m/s

load N_direction1.mat
distance1=distance_TCA;
load T_direction1.mat
distance2=distance_TCA;
load W_direction1.mat
distance3=distance_TCA;
load times.mat

plot(time,distance1)
hold on
plot(time,distance2)
hold on
plot(time,distance3)
xlabel('Deflection times (day)')
ylabel('Distance from Earth at TCA (Re)')
legend('N','T','W')
title('Deflection with 1 m/s')