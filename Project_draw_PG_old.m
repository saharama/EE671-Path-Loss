close all;
clear;

figure; hold on; %grid on;

plot([0, 500], [70, 70], '-o')
plot([20, 20], [0, 500], '-o')
plot(0,0, 'o')
plot([50, 50, 600], [600, 110, 110], '-o')
plot([50, 50, 600], [-10, 50, 50], '-o')

%1 case, 1 point, only direct
%path gain is equal to (wavelength/4pi)^2 ((1/r)e^-jbetar)^2
%because transmitter is 0, 0, makes calculating easier

%todo: put path gain in own function?
%      how does diffraction work again lol

wavelength = 3*10^8/1000000000; %c/f
beta = 2*pi/wavelength;
%x = 0; y = 70; R = sqrt(x^2 + y^2);
%PG = (wavelength/4/pi)^2*abs(exp(-1i*beta*R)/R)^2;
n = 1;
Q1x = 50; Q1y = 50; Q2x = 50; Q2y = 110;

ab = [0:.05:500;70*ones(1,10001)];
Rfull = sqrt(ab(1,:).^2 + ab(2,:).^2);
incidfull = exp(-1i*beta.*Rfull)./Rfull;

PGfull = (wavelength/4/pi)^2*abs(incidfull).^2;

%AB, to 70 since Q1 blocks at that point
for x = 0:.05:500
    y = 70;
    xAB(n) = x;
    if x < 70
        R = sqrt(x^2 + y^2);
        incid = exp(-1i*beta*R)/R;
%         PG = (wavelength/4/pi)^2*abs(exp(-1i*beta*R)/R)^2;
%         PGAB(n) = PG;
    else
        incid = 0;
    end
    incidAB(n) = incid;
    
    %diff1 - diffraction due to Q1
%     sp = sqrt(50^2 + 50^2); s = sqrt((Q1x - x)^2 + (Q1y - y)^2); phip = 90; 
%     if x < 70
%         phi = 180 - acos(20/s)/2/pi*180;
%     elseif x == 70
%         phi = 180;
%     else
%         phi = 180 + acos(20/s)/2/pi*180;
%     end
%     L = 1/wavelength*sp*s/(s+sp);
%     [ds,dh] = wdc(L, phi, phip, 90, 1.5);
%     diff1 = -sqrt(wavelength)*ds/sqrt(sp*s*(sp+s))*exp(-1i*beta*(sp+s));
%     diff1AB(n) = diff1;
    
    %diff2 - diffraction due to Q2
    %sp = sqrt(50^2+110^2); s = sqrt((Q2x - x)^2 + (Q2y - y)^2);
    
    %cw
    %angle1 = 90 + asin(50/sp)/2/pi*180; angle2 = 90 + asin(50/s)/2/pi*180;
    
    %ccw
    %angle1 = 
    
    PG = (wavelength/4/pi)^2*abs(incid)^2;
    PGAB(n) = PG;
    n = n+1;
end

%CD
n = 1;
for y = 0:.05:500
    yCD(n) = y;
    x = 20;
    R = sqrt(x^2 + y^2);
    PG = (wavelength/4/pi)^2*abs(exp(-1i*beta*R)/R)^2;
    PGCD(n) = PG;
    
    n = n+1;
end

figure
PGABdB = 10*log10(PGAB);
plot(xAB, PGAB);

figure
PGfulldB = 10*log10(PGfull);
plot(ab(1,:), PGfull);

figure
PGCDdB = 10*log10(PGCD);
plot(yCD, PGCDdB);