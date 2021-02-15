function [a, A] = Elevation_and_Azimuth(r, R, lstOS, Lat)
%{
  Elevation_and_Azimuth function calculates the elevation angle and azimuth
  of the satellite in the topocentric horizon coordinate system located at
  the observation site.

    a     - Elevation angle of the satellite in radians
    A     - Azimuth of the satellite in radians
    r     - Position vector of the satellite in geocentric equatorial frame
    R     - Position vector of the observation site in geocentric 
            equatorial frame
    lstOS - Sidereal time of the observation site in degrees
    Lat   - Latitude of the observation site in degrees
    
  User M-function required : None
  User subfunction required: None

  Author: Tahsin Çaðrý Þiþman
%}

% The relative position vector of the satellite with respect to the
% observation site in the geocentric equatorial frame:
rho = r - R;

% Transformation matrix from the geocentric equatorial frame to the
% topocentric horizon frame:
th = deg2rad(lstOS);
ph = deg2rad(Lat);

Q_Xx = [-sin(th)           cos(th)          0 
        -sin(ph)*cos(th)  -sin(ph)*sin(th)  cos(ph)
         cos(ph)*cos(th)   cos(ph)*sin(th)  sin(ph)];
    
% Writing the relative position vector in the topocentric horizon frame:
    rhoTH = (Q_Xx*rho')';

% The unit vector in the direction of the relative position vector:
    rhoTHdir = rhoTH/norm(rho);

% Finding the elevation for the satellite in topocentric horizon frame
% located at the observation site:
    a = asin(rhoTHdir(3));

% Finding the azimuth for the satellite in topocentric horizon frame
% located at the observation site:
    A = acos(rhoTHdir(2)/cos(a));
    if sign(sin(A)) == sign(rhoTHdir(1)/cos(a))
    
    else
        A = 2*pi - A;
    end
    
end % Elevation_and_Azimuth