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
Q1x = 50; Q1y = 50; Q2x = 50; Q2y = 110;

%AB direct field
ab = [0:.05:500;70*ones(1,10001)];
RAB = sqrt(ab(1,:).^2 + ab(2,:).^2);
incidAB = exp(-1i*beta.*RAB)./RAB;
incidAB(1001:end) = 0; %50/.05 = 1000

%diff1 AB- diffraction due to Q1
sp = sqrt(50^2 + 50^2);
s = sqrt((Q1x - ab(1,:)).^2+(Q1y - ab(2,:)).^2);
L = 1/wavelength.*sp.*s./(s+sp);
v1 = [zeros(1,10001);-1*ones(1,10001)]; phip = 45;
v2 = [ab(1,:) - 50; ab(2,:) - 50]./s;
phi = acos(dot(v1,v2))*180/pi;
phi(1001:end) = 360 - phi(1001:end);

for n = 1:1:length(phi)
    [ds,dh] = wdc(L(n),phi(n),phip,90,1.5);
    dsAB(n) = ds;
end

diff1AB = -sqrt(wavelength).*dsAB./sqrt(sp.*s.*(sp+s)).*exp(-1i*beta.*(sp+s));

%diff2 AB - diffraction due to Q2
sp = sqrt(50^2 + 110^2);
s = sqrt((Q2x - ab(1,:)).^2+(Q2y - ab(2,:)).^2);
L = 1/wavelength.*sp.*s./(s+sp);
v1 = [zeros(1,10001);ones(1,10001)];
v2 = [ab(1,:) - 50; ab(2,:) - 110]./s;
phip = acos(dot([50,110]/sqrt(110^2+50^2),[0,1]))*180/pi;
phi = acos(dot(v1,v2))*180/pi;
phi(1001:end) = 360 - phi(1001:end);

for n = 1:1:length(phi)
    [ds,dh] = wdc(L(n),phi(n),phip,90,1.5);
    dsAB(n) = ds;
end
diff2AB = -sqrt(wavelength).*dsAB./sqrt(sp.*s.*(sp+s)).*exp(-1i*beta.*(sp+s));

PGAB = (wavelength/4/pi)^2 .* abs(incidAB + diff1AB + diff2AB).^2;

%CD
cd = [20*ones(1,10001);0:.05:500];

%CD direct field
RCD = sqrt(cd(1,:).^2 + cd(2,:).^2);
incidCD = exp(-1i*beta.*RCD)./RCD;

%diff1 CD- diffraction due to Q1
sp = sqrt(50^2 + 50^2);
s = sqrt((Q1x - cd(1,:)).^2+(Q1y - cd(2,:)).^2);
L = 1/wavelength.*sp.*s./(s+sp);
v1 = [zeros(1,10001);-1*ones(1,10001)]; phip = 45;
v2 = [cd(1,:) - 50; cd(2,:) - 50]./s;
phi = acos(dot(v1,v2))*180/pi;
phi(1001:end) = 360 - phi(1001:end);

for n = 1:1:length(phi)
    [ds,dh] = wdc(L(n),phi(n),phip,90,1.5);
    dsCD(n) = ds;
end

diff1CD = -sqrt(wavelength).*dsCD./sqrt(sp.*s.*(sp+s)).*exp(-1i*beta.*(sp+s));

%diff2 CD - diffraction due to Q2
sp = sqrt(50^2 + 110^2);
s = sqrt((Q2x - cd(1,:)).^2+(Q2y - cd(2,:)).^2);
L = 1/wavelength.*sp.*s./(s+sp);
v1 = [zeros(1,10001);ones(1,10001)];
v2 = [cd(1,:) - 50; cd(2,:) - 110]./s;
phip = acos(dot([50,110]/sqrt(110^2+50^2),[0,1]))*180/pi;
phi = acos(dot(v1,v2))*180/pi;
phi(1001:end) = 360 - phi(1001:end);

for n = 1:1:length(phi)
    [ds,dh] = wdc(L(n),phi(n),phip,90,1.5);
    dsCD(n) = ds;
end
diff2CD = -sqrt(wavelength).*dsCD./sqrt(sp.*s.*(sp+s)).*exp(-1i*beta.*(sp+s));

PGCD = (wavelength/4/pi)^2 .* abs(incidCD + diff1CD + diff2CD).^2;

figure
PGABdB = 10*log10(PGAB);
plot(ab(1,:), PGABdB);
title('Path Gain on a Street')
xlabel('x(m)');
ylabel('Path Gain (dB)');

figure
PGCDdB = 10*log10(PGCD);
plot(cd(2,:), PGCDdB);
title('Path Gain on a Street')
xlabel('y(m)');
ylabel('Path Gain (dB)');
