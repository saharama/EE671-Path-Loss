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

wavelength = 1/1000000000;
beta = 2*pi/wavelength;
x = 0; y = 70; d = sqrt(x^2 + y^2);
PG = (wavelength/4/pi)^2*abs(exp(-1i*beta*d)/d)^2;
n = 1;

%AB, to 70 since Q1 blocks at that point
for x = 0:.05:70
    xAB(n) = x;
    y = 70;
    d = sqrt(x^2 + y^2);
    PG = (wavelength/4/pi)^2*abs(exp(-1i*beta*d)/d)^2;
    PGAB(n) = PG;
    
    n = n+1;
end

%CD
n = 1;
for y = 0:.05:500
    yCD(n) = y;
    x = 20;
    d = sqrt(x^2 + y^2);
    PG = (wavelength/4/pi)^2*abs(exp(-1i*beta*d)/d)^2;
    PGCD(n) = PG;
    
    n = n+1;
end

figure
PGABdB = 10*log10(PGAB);
plot(xAB, PGABdB);

figure
PGCDdB = 10*log10(PGCD);
plot(yCD, PGCDdB);