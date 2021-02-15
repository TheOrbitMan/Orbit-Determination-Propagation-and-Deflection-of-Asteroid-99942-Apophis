%-------------------------------------------------------------------------%
% Computes right ascension and declination from the state vector
%
% REFERENCE:
% Curtis, H. D.(2013) Orbital Mechanics for Engineering Students, 3. Edition, Elsevier.
%
% Author: Þahin Ulaþ KÖPRÜCÜ 
%-------------------------------------------------------------------------%
function [RA,DEC] = RA_DEC_from_sv(r)

rmag=norm(r);

l=r(1)/rmag;m=r(2)/rmag;n=r(3)/rmag;

DEC=asind(n);

if m>0
    RA=acosd(l/cosd(DEC));
else
    RA=360-acosd(l/cosd(DEC));
end


end

