%-------------------------------------------------------------------------%
% Computes two-body acceleration
% 
% Author: Þahin Ulaþ KÖPRÜCÜ
%-------------------------------------------------------------------------%
function [rate] = kepler(pos,mu)

accel=(-mu*pos)/norm(pos)^3;
rate=[accel(1) accel(2) accel(3)];

end

