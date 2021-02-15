%-------------------------------------------------------------------------%
% Computes two-body acceleration
% 
% Author: �ahin Ula� K�PR�C�
%-------------------------------------------------------------------------%
function [rate] = kepler(pos,mu)

accel=(-mu*pos)/norm(pos)^3;
rate=[accel(1) accel(2) accel(3)];

end

