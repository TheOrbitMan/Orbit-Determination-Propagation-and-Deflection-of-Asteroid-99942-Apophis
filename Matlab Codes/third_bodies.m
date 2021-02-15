%-------------------------------------------------------------------------%
% Computes acceleration due to Third-Bodies
% 
% REFERENCE:
% Curtis, H. D.(2013) Orbital Mechanics for Engineering Students, 3. Edition, Elsevier.
%
% Author: Þahin Ulaþ KÖPRÜCÜ
%-------------------------------------------------------------------------%
function [rate] = third_bodies(pos,JD)

%AU(km)
au=149597870.700;
%au=149597871;

%Heliocentric position vector of third-bodies.(km)
[r_mercury, v_mercury] = planet_elements_and_sv(1,JD);
[r_venus, v_venus] = planet_elements_and_sv(2,JD);
[r_earth, v_earth] = planet_elements_and_sv(3,JD);
[r_mars, v_mars] = planet_elements_and_sv(4,JD);
[r_jupiter, v_jupiter] = planet_elements_and_sv(5,JD);
moon_req = lunar_position(JD);
moon_veq = [0 0 0];
[moon_rec,moon_vec] = Eq2Ec(moon_req,moon_veq,JD);

%Heliocentric position vector of third-bodies.(AU)
R_mercury=r_mercury/au;
R_venus=r_venus/au;
R_earth=r_earth/au;
R_mars=r_mars/au;
R_jupiter=r_jupiter/au;
r_moon=moon_rec+r_earth;
R_moon=r_moon/au;

%Gravitational parameters for the planets(From GMAT)
mu_mercury=22032.080486418*(86400^2)/(au^3);
mu_venus=324858.59882646*(86400^2)/(au^3);
mu_earth=398600.4415*(86400^2)/(au^3);
mu_mars=42828.314258067*(86400^2)/(au^3);
mu_jupiter=126712767.8578*(86400^2)/(au^3);
mu_moon=4902.8005821478*(86400^2)/(au^3);

%Third-bodies accelerations.(au/day^2)
accel1=mu_mercury*((R_mercury-pos)/norm(R_mercury-pos)^3-R_mercury/norm(R_mercury)^3);
accel2=mu_venus*((R_venus-pos)/norm(R_venus-pos)^3-R_venus/norm(R_venus)^3);
accel3=mu_earth*((R_earth-pos)/norm(R_earth-pos)^3-R_earth/norm(R_earth)^3);
accel4=mu_mars*((R_mars-pos)/norm(R_mars-pos)^3-R_mars/norm(R_mars)^3);
accel5=mu_jupiter*((R_jupiter-pos)/norm(R_jupiter-pos)^3-R_jupiter/norm(R_jupiter)^3);
accel6=mu_moon*((R_moon-pos)/norm(R_moon-pos)^3-R_moon/norm(R_moon)^3);

accel=accel1+accel2+accel3+accel4+accel5+accel6;

rate=[accel(1) accel(2) accel(3)];
end

