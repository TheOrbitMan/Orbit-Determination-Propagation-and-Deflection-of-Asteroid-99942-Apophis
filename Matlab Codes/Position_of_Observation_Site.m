function [R, lstOS] = Position_of_Observation_Site(date,Observation_Site)
%{
  Position_of_Observation_Site function calculates the position vector of
  the observation site in geocentric equatorial frame.

    R                - Position of the observation site in geocentric
                       equatorial frame
    Observation_Site - Involves east longitude, latitude and altitude of
                       the observation site
    date             - Date of the instant
    UT               - Universal time of the instant
    EL               - East longitude of the observation site in degrees
    Lat              - Latitude of the observation site in degrees
    lstOS            - Sidereal time of the observation site in degrees
    H                - Altitude of the observation site in km
    
  User M-function required : LST.m (of Curtis)
  User subfunction required: None

  Author: Tahsin Çaðrý Þiþman

%}

% The radius of Earth RE in km units and the oblateness of Earth which is
% dimensionless are:
RE = 6378;
f = 0.00335;

% Observation site:
EL = Observation_Site(1);
Lat = Observation_Site(2);
H = Observation_Site(3);

% Universal time of the instant:
UT = date(4) + date(5)/60 + date(6)/3600;

% Local side real time for the observation site (determining the angular
% position of the observation site in the geocentric equatorial frame). 
% Note that the local sidereal time output of the LST function is in the 
% unit degree.
lstOS = LST(date(1),date(2),date(3),UT,EL);

% Finding the position vector for the observation site in the geocentric
% equatorial frame by using (5.56) of Curtis, 3rd Ed. 
th = deg2rad(lstOS);
ph = deg2rad(Lat);

Rx = (RE/sqrt(1-(2*f-f^2)*sin(ph)^2) + H)*cos(ph)*cos(th);
Ry = (RE/sqrt(1-(2*f-f^2)*sin(ph)^2) + H)*cos(ph)*sin(th);
Rz = (RE*(1-f)^2/sqrt(1-(2*f-f^2)*sin(ph)^2) + H)*sin(ph);

R = [Rx, Ry, Rz];
    
end % Position_of_Observation_Site


