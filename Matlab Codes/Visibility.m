%-------------------------------------------------------------------------%
% Calculation of naked-eye visibility of the Apophis from a location on Earth
%
% Author: Þahin Ulaþ KÖPRÜCÜ 
%-------------------------------------------------------------------------%

clc;clear;

%AU
au=149597871; %km

%Observation site (ANKARA)
EL=32.8597;   %(East longitude) degree
LAT=39.9334;  %degree
H=0.9;        %km
Observation_Site=[EL,LAT,H];

%Time of closest approach(TCA) (2029-Apr-13 21:46:12.6856)
TCA =2462240.407091269; %JD
date=[2029 04 13 21 46 12.6856];

%Position vector of the Asteroid in Sun Centered Ecliptic frame at TCA (au)
%Precise ephemeris: https://ssd.jpl.nasa.gov/horizons.cgi#results
R_asteroid=[-9.175061798507436E-01 -4.050876512426926E-01 6.998818320904296E-05];

%Position vector of the Earth in Sun Centered Ecliptic frame at TCA (au)
%Precise ephemeris: https://ssd.jpl.nasa.gov/horizons.cgi#results
R_earth=[-9.173784794862199E-01 -4.053035544141292E-01 2.967146247411548E-05];

%Position vector of the Earth in Sun Centered Ecliptic frame at TCA (au)
%R_earth=planet_elements_and_sv(3,TCA)/au;

%Position vector of Asteroid relative to Earth in Sun Centered Ecliptic frame(au)
r_rel_ec=R_asteroid-R_earth;

%Asteroid-Earth distance at TCA (au)
distace_TCA=norm(r_rel_ec);

%Position vector of Asteroid relative to Earth in Equatorial frame(au)
[r_rel_eq,v_rel_eq] = Ec2Eq(r_rel_ec,[0 0 0],TCA);

%Position vector of Sun relative to Earth in Equatorial frame(au)
[r_sun_eq,v_sun_eq] = Ec2Eq(-R_earth,[0 0 0],TCA);

%Position vector of Observer in Equatorial frame(km)
[r_site_eq, lstOS] = Position_of_Observation_Site(date,Observation_Site);

%Azimuth and Elevation of Asteroid
[a, A] = Elevation_and_Azimuth(r_rel_eq*au, r_site_eq, lstOS, LAT);
el_asteroid=rad2deg(a);
az_asteroid=rad2deg(A);

%Azimuth and Elevation of Sun
[a, A] = Elevation_and_Azimuth(r_sun_eq*au, r_site_eq, lstOS, LAT);
el_sun=rad2deg(a);
az_sun=rad2deg(A);

%Checking for Shadow condition (Cylindrical shadow model)
light_switch = los(r_rel_eq*au, r_sun_eq*au);

%Appropriate conditions for the naked eye observation
if el_asteroid>10 && el_sun<-6 && light_switch==1
    fprintf('Visible');
    el_asteroid
    az_asteroid
    el_sun
else
    fprintf('Not Visible');
end