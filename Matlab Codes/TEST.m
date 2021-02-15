%-------------------------------------------------------------------------%
% Test.m 
%
% This code is written within the scope of the conference paper: 
% Koprucu, S. U. (2019) Orbit Determination, Propagation and Deflection of
% Asteroid (99942) Apophis, 10th Ankara International Aerospace Conference.
%
% Author: Þahin Ulaþ KÖPRÜCÜ 
%-------------------------------------------------------------------------%
clc;clear;
format long

%% PRELIMINARY ORBIT DETERMINATION (POD)

% Observations are taken from:
% https://minorplanetcenter.net/db_search/show_object?utf8=%E2%9C%93&object_id=apophis------>download------>https://minorplanetcenter.net/tmp/apophis.txt
% Please enter the OBSERVATIONS.txt file with proper format.
%-------------------------------------------------------------------------%
obs =fopen('OBSERVATIONS.txt');
line=fgetl(obs);
i=1;
while ischar(line)
    year(i) =str2num(line(16:19));
    month(i)=str2num(line(21:22));
    day(i)  =str2num(line(23:25));
    dayf    =str2num(line(26:32));
    hf      =dayf*24;
    mf      =(hf-floor(hf))*60;
    sf      =(mf-floor(mf))*60;
    hour(i) =floor(hf);
    minute(i)=floor(mf);
    second(i)=sf;
    jd(i)   =J0(year(i),month(i),day(i))+(hour(i)+minute(i)/60+second(i)/3600)/24;
    rhour   =str2num(line(33:35));
    rmin    =str2num(line(36:38));
    rsec    =str2num(line(39:44));
    RA(i)   =((rhour*3600+rmin*60+rsec)*15)/3600;
    sign    =(line(45:45));
    ddeg    =str2num(line(46:47));
    dmin    =str2num(line(48:50));
    dsec    =str2num(line(51:56));
    if sign=='-'
    DEC(i)  =(ddeg+(dmin)/60+(dsec)/3600)*(-1);
    else
    DEC(i)  =(ddeg+(dmin)/60+(dsec)/3600);
    end
    line    = fgetl(obs);
    i=i+1;
end
fclose(obs);
%-------------------------------------------------------------------------%

% Constants
au=149597871;                        % km
mu=((1.327*10^11)*(86400^2)/(au^3)); % au^3/day^2 (Gravitational parameter for Sun)
Re=6378;                             % km (Radius of Earth)
j2000=2400000.5;                     % day

%1) Preliminary Orbit Determination(Gauss and Laplace)
%-------------------------------------------------------------------------%
% Sun to Earth vector in ecliptic frame with 3 different observation epoch.(km and km/s)
[r1, v1] = planet_elements_and_sv(3,jd(1));
[r2, v2] = planet_elements_and_sv(3,jd(2));
[r3, v3] = planet_elements_and_sv(3,jd(3));

% Ecliptic frame to equatorial frame.(km and km/s)
[R1,V1] = Ec2Eq(r1,v1,jd(1));
[R2,V2] = Ec2Eq(r2,v2,jd(2));
[R3,V3] = Ec2Eq(r3,v3,jd(3));

% Sun to Earth vector in equatorial frame.(au)
R1=R1/au;R2=R2/au;R3=R3/au;

% Earth to Asteroid line of sight unit vectors.
rhohat1=[cosd(DEC(1))*cosd(RA(1)) cosd(DEC(1))*sind(RA(1)) sind(DEC(1))];
rhohat2=[cosd(DEC(2))*cosd(RA(2)) cosd(DEC(2))*sind(RA(2)) sind(DEC(2))];
rhohat3=[cosd(DEC(3))*cosd(RA(3)) cosd(DEC(3))*sind(RA(3)) sind(DEC(3))];

% Preliminary Orbit Determination (GAUSS's and LAPLACE's methods)
xi=2; %Initial guess for r2.(au)
% State vector at second observation epoch in the equatorial frame from LAPLACE's method.(au,au/day)
[r2LAPLACE, v2LAPLACE] = LAPLACE(rhohat1, rhohat2, rhohat3, R1, R2, R3, jd(1), jd(2), jd(3), mu, xi);
% State vector at second observation epoch in the equatorial frame from GAUSS's method.(au,au/day)
[r2GAUSS, v2GAUSS] = GAUSS(rhohat1, rhohat2, rhohat3, R1, R2, R3, jd(1), jd(2), jd(3), mu, xi);

% Orbital elements in the ecliptic frame
[rec1,vec1] = Eq2Ec(r2GAUSS,v2GAUSS,jd(2));
coe1        = coe_from_sv(rec1,vec1,mu);
e1=coe1(1)
raan1=radtodeg(coe1(2))
i1=radtodeg(coe1(3))
w1=radtodeg(coe1(4))
ta1=radtodeg(coe1(5))
a1=coe1(6)
%-------------------------------------------------------------------------%

%% ORBIT ESTIMATION
%-------------------------------------------------------------------------%
%Nominal vectors obtained from POD in equatorial frame.(r2GAUSS-v2GAUSS or r2LAPLACE-v2LAPLACE)
xnom(1,1) = r2LAPLACE(1);
xnom(2,1) = r2LAPLACE(2);
xnom(3,1) = r2LAPLACE(3);
xnom(4,1) = v2LAPLACE(1);
xnom(5,1) = v2LAPLACE(2);
xnom(6,1) = v2LAPLACE(3);
%Iteration of the Differential Correction method.
noiseRA=0;  %RA observation noise in deg.
noiseDEC=0; %DEC observation noise in deg.
secondobs=2;lastobs=i-1;
criteria=10e-2;rel_error=100;
rel=zeros(6,1);
while (rel_error)>criteria
percentchg = 0.0001;
[atwa, atwb]=findatwaatwb(secondobs, lastobs, RA, DEC, jd, year, month, day, hour, minute, second, percentchg, xnom, mu, au, noiseRA, noiseDEC);
P=inv(atwa);
dx=P*atwb;
rel_error=max(abs((xnom+dx)-xnom)./abs(xnom+dx))*100;
xnom=xnom+dx;
end
%Improved state vectors and orbital elements.
rimpr=[xnom(1) xnom(2) xnom(3)]';vimpr=[xnom(4) xnom(5) xnom(6)]';
[rec2,vec2] = Eq2Ec(rimpr,vimpr,jd(2));
coe2 = coe_from_sv(rec2,vec2,mu);
e2=coe2(1)
raan2=radtodeg(coe2(2))
i2=radtodeg(coe2(3))
w2=radtodeg(coe2(4))
ta2=radtodeg(coe2(5))
a2=coe2(6)
%-------------------------------------------------------------------------%

%% ORBIT PROPAGATION 
%-------------------------------------------------------------------------%

%For use in more precise propagation with state vector from NASA JPL HORIZONS
au=149597870.700;
mu=132712440017.99*((86400^2)/(au^3));

%Time of closest approach(julian day)
%tca=2.462240400000000e+06; %ESA (http://neo.ssa.esa.int/search-for-asteroids?tab=closeapp&des=99942%20Apophis)
tca=2462240.407091269;     %NASA (https://ssd.jpl.nasa.gov/sbdb.cgi?sstr=99942;orb=0;cov=0;log=0;cad=1#cad)	

%Propagation duration(day)
tf=tca-jd(2);

%Initial state determined from the orbit determination at second observation epoch
state=[rec2 vec2];

%Initial state from NASA JPL HORIZONS in Sun centered MJ2000Ec at second observation epoch (https://ssd.jpl.nasa.gov/horizons.cgi#results)
%state=[-4.648982993233945e-01 9.443590464045293e-01 -6.125110499611269e-02 -1.467312309906414e-02 -5.015572583174898e-03 -8.744338994778924e-05];
%state=[-4.648982685665431E-01 9.443590565730078E-01 -6.125110430798041E-02 -1.467312332529396E-02 -5.015572021998710E-03 -8.744361686909425E-05]; %Old version

%FIXED STEPSIZE 
%Time step(day)
dt=0.5;JD=jd(2);
i=1;
while JD<tca
    [state] = RK4(state,min(dt,tca-JD),mu,JD);
%    [state] = RK8(state,min(dt,tca-JD),mu,JD);
%    [state] = PrinceDormand(state,min(dt,tca-JD),mu,JD);
    RX(i)=state(1);RY(i)=state(2);RZ(i)=state(3);
    VX(i)=state(4);VY(i)=state(5);VZ(i)=state(6);
    R=[state(1) state(2) state(3)];
    JD=JD+min(dt,tca-JD);
    %Position of Earth wrt Sun in ecliptic frame(km)
    R_earth=planet_elements_and_sv(3,JD);
    %Distance between Asteroid-Earth(au)
    dist(i)=norm(R-R_earth/au);
    t(i)=dt*i;
    JD_array(i)=JD;
    i=i+1;
end
plot(t./365,dist)
xlabel('Time since epoch (year)')
ylabel('Distance from Earth (AU)')

r_final=[state(1) state(2) state(3)];
[re,ve]=planet_elements_and_sv(3,tca);
distance_at_TCA=norm(r_final-(re/au))

%ADAPTIVE STEPSIZE ODE45 TO CHECK
%{
tspan=[jd(2) tca];
f0=state';
options=odeset('RelTol',1e-13,'MaxStep',1);
[t,f]= ode45(@(t,y) rates(t,y,mu), tspan, f0, options);
RX=f(:,1);RY=f(:,2);RZ=f(:,3);
VX=f(:,4);VY=f(:,5);VZ=f(:,6);
r_final=[f(length(f),1) f(length(f),2) f(length(f),3)];
[r_earth,v_earth]=planet_elements_and_sv(3,tca); %km
distance_at_TCA=norm(r_final-(r_earth/au))       %au
%}

%-------------------------------------------------------------------------%

%% ORBIT DEFLECTION
%For the orbit deflection, it is recommended to use initial state vector
%from NASA JPL HORIZONS because it will give closer result to reference values
%-------------------------------------------------------------------------%
%Time interval between the deflections(day)
t_defl=5;
%Number of points
fact=100;
%Radius of Earth(km)
Re_radius=6378;
%Magnitude of the velocity change of asteroid
mag_delv=1;                    %m/s
ms2auday=86400/(au*1000);
delv=mag_delv*ms2auday;        %au/day
%Direction of the velocity change due to deflection (It is defined in NTW frame)
delv=[0 delv 0];

%State vector of the Earth wrt Sun in ecliptic frame at tca(km km/s)
[re,ve] = planet_elements_and_sv(3,tca);

if mod(t_defl,dt)~=0
    error('t_defl has to be divided to timestep exactly');
end

%Please choose the same propagator as in the previous part

for i=1:fact
     [Q] = SCI2NTW(RX(length(t)-(t_defl/dt)*i),RY(length(t)-(t_defl/dt)*i),RZ(length(t)-(t_defl/dt)*i),...
         VX(length(t)-(t_defl/dt)*i),VY(length(t)-(t_defl/dt)*i),VZ(length(t)-(t_defl/dt)*i));
     delv_t=Q'*delv';
     R=[RX(length(t)-(t_defl/dt)*i) RY(length(t)-(t_defl/dt)*i) RZ(length(t)-(t_defl/dt)*i)];
     V=[VX(length(t)-(t_defl/dt)*i) VY(length(t)-(t_defl/dt)*i) VZ(length(t)-(t_defl/dt)*i)]+delv_t';
     state=[R V];
     JD=JD_array(length(t)-(t_defl/dt)*i);
     while JD<tca
    [state] = RK4(state,min(dt,tca-JD),mu,JD);
%    [state] = RK8(state,min(dt,tca-JD),mu,JD);
%    [state] = PrinceDormand(state,min(dt,tca-JD),mu,JD);
     JD=JD+min(dt,tca-JD);
     end   
     R=[state(1) state(2) state(3)];
     distance_TCA(i)=((norm(R-re/au))*au)/Re_radius;
     time(i)=tca-JD_array(length(t)-(t_defl/dt)*i);
end
%  save('T_direction.mat','distance_TCA','-double');
%  save('times.mat','time','-double');
plot(distance_TCA)

%Plot of the Apophis-Sun distance
%{
for i=1:fact
R=[RX(length(t)-(t_defl/dt)*i) RY(length(t)-(t_defl/dt)*i) RZ(length(t)-(t_defl/dt)*i)];
r_ast(i)=norm(R);
time_left(i)=tca-JD_array(length(t)-(t_defl/dt)*i);
end
plot(time_left,r_ast)
xlabel('Deflection times (day)')
ylabel('Distance from Sun (AU)')
%}