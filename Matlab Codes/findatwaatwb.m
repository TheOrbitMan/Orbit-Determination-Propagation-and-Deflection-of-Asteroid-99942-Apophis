%-------------------------------------------------------------------------%
% Differential Correction
%
% REFERENCE:
% Vallado, D. A.(2013) Fundamentals of Astrodynamics and Applications, 4. Edition, Microcosm Press.
%
% This code is taken from:
% https://celestrak.com/software/vallado-sw.php-------->Computer software in Matlab
% But it is modified.
%
% Last modified: Þahin Ulaþ KÖPRÜCÜ
%-------------------------------------------------------------------------%
function [atwa, atwb] = findatwaatwb(secondobs, lastobs, RA, DEC, jd, year, month, day, hour, minute, second, percentchg, xnom, mu, au, noiseRA, noiseDEC)
    
    rad=pi/180;
    
    sumRA  = 0;
    sumDEC = 0;
    
    % ------------- reset these since they will accumulate ---------
    atwa=zeros(6,6);
    atwb=zeros(6,1);
    % ------------------- loop through all the observations ------------------
    for i = secondobs: (lastobs-1)
        
        % --------- propagate the nominal vector to the epoch time -----------
        dt=jd(i+1)-jd(2); %day
        rnom(1) = xnom(1,1);
        rnom(2) = xnom(2,1);
        rnom(3) = xnom(3,1);
        vnom(1) = xnom(4,1);
        vnom(2) = xnom(5,1);
        vnom(3) = xnom(6,1);
        [rec,vec] = Eq2Ec(rnom,vnom,jd(2));
        [r,v] = Two_Body_Propagator(rec,vec,dt,mu);
        [req,veq] = Ec2Eq(r,v,jd(i+1));
        
        % ------------------------- find b matrix ----------------------------
        %Computed RA and DEC
        [R] = solar_position(year(i+1),month(i+1),day(i+1),hour(i+1),minute(i+1),second(i+1));
        Req=-R'/au;
        [compRA,compDEC] = RA_DEC_from_sv(req-Req);
        
        %Observed RA and DEC
        obsRA=RA(i+1);
        obsDEC=DEC(i+1);
        
        b(1,1)=(obsRA-compRA)*rad;
        b(2,1)=(obsDEC-compDEC)*rad;
     
        % ------------------------ find a matrix ----------------------------
           
        % ----- perturb each element in the state (elements or vectors) ------
        for j= 1:6
            [deltaamt, xnomp] = finitediff(j, percentchg, xnom);
            rnomp(1) = xnomp(1,1);
            rnomp(2) = xnomp(2,1);
            rnomp(3) = xnomp(3,1);
            vnomp(1) = xnomp(4,1);
            vnomp(2) = xnomp(5,1);
            vnomp(3) = xnomp(6,1);
            [recp,vecp] = Eq2Ec(rnomp,vnomp,jd(2));
            [rp,vp] = Two_Body_Propagator(recp,vecp,dt,mu);
            [reqp,veqp] = Ec2Eq(rp,vp,jd(i+1));
            [compRAp,compDECp] = RA_DEC_from_sv(reqp-Req);
            a(1,j)=((compRAp-compRA)*rad)/(deltaamt);
            a(2,j)=((compDECp-compDEC)*rad)/(deltaamt);
        end
        
        at = a';
        
        % ------------------------- assign weights ---------------------------
        if noiseRA==0 && noiseDEC==0
            atw=at;
        else
        w1=1/(noiseRA*rad)^2;
        w2=1/(noiseDEC*rad)^2;
        w=[w1 0;0 w2];
        atw=at*w;
        end
        
        % ----------------- find the atwa / atwb matrices --------------------
        atwaacc = atw * a;
        atwbacc = atw * b;
        % ------------------- accumulate the matricies -----------------------
        atwa=atwa+atwaacc;
        atwb=atwb+atwbacc; 
        
    end