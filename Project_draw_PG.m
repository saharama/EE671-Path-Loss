close all; clear;
figure; hold on; grid on;

%plot street intersection
plot([0, 500], [70, 70], '-o')
plot([20, 20], [0, 500], '-o')
plot(0,0, 'o')
plot([50, 50, 600], [600, 110, 110], '-o')
plot([50, 50, 600], [-10, 50, 50], '-o')
xlabel('x(m)'); ylabel('y(m)');
title('Street');

%whole parameters
wavelength = 3*10^8/1000000000; %c/f
beta = 2*pi/wavelength;
Q1x = 50; Q1y = 50; Q2x = 50; Q2y = 110;

%AB direct field
ab = [0:.05:500;70*ones(1,10001)];
RAB = sqrt(ab(1,:).^2 + ab(2,:).^2);
incidAB = exp(-1i*beta.*RAB)./RAB;
incidAB(1401:end) = 0; %50/.05 = 1000

%diff1 AB- diffraction due to Q1
sp = sqrt(50^2 + 50^2); %calculate distance from transmitter to corner
s = sqrt((Q1x - ab(1,:)).^2+(Q1y - ab(2,:)).^2); %calculate distance from corner to receiver
L = 1/wavelength.*sp.*s./(s+sp); %calculate L
v1 = [zeros(1,10001);-1*ones(1,10001)]; phip = 45; %define wall unit vector, phi'
v2 = [ab(1,:) - 50; ab(2,:) - 50]./s; %calculate unit vectors from corner to receiver
phi = acos(dot(v1,v2))*180/pi; %find angle
phi(1001:end) = 360 - phi(1001:end); %change negative (=pos) cosines to positives

%take diffraction constant from function
for n = 1:1:length(phi)
    [ds,dh] = wdc(L(n),phi(n),phip,90,1.5);
    dsAB(n) = ds;
end

%calculate path gain based on summed matrix field
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

%take diffraction constant from function
for n = 1:1:length(phi)
    [ds,dh] = wdc(L(n),phi(n),phip,90,1.5);
    dsAB(n) = ds;
end
diff2AB = -sqrt(wavelength).*dsAB./sqrt(sp.*s.*(sp+s)).*exp(-1i*beta.*(sp+s));

PGAB = (wavelength/4/pi)^2 .* abs(incidAB + diff1AB + diff2AB).^2;
PGAB_d = (wavelength/4/pi)^2 .* abs(diff1AB + diff2AB).^2;
PGAB_i = (wavelength/4/pi)^2 .* abs(incidAB).^2;

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

%take diffraction constant from function
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

%take diffraction constant from function
for n = 1:1:length(phi)
    [ds,dh] = wdc(L(n),phi(n),phip,90,1.5);
    dsCD(n) = ds;
end
diff2CD = -sqrt(wavelength).*dsCD./sqrt(sp.*s.*(sp+s)).*exp(-1i*beta.*(sp+s));

%calculate path gain from summed matrices
PGCD = (wavelength/4/pi)^2 .* abs(incidCD + diff1CD + diff2CD).^2;
PGCD_i = (wavelength/4/pi)^2 .* abs(incidCD).^2;
PGCD_d = (wavelength/4/pi)^2 .* abs(diff1CD + diff2CD).^2;

%plot path gain from east to west
figure; hold on;
PGABdB = 10*log10(PGAB);
PGAB_idB = 10*log10(PGAB_i);
PGAB_ddB = 10*log10(PGAB_d);
plot(ab(1,:), PGABdB,'DisplayName','total');
plot(ab(1,:), PGAB_idB,'DisplayName','incident');
plot(ab(1,:), PGAB_ddB,'DisplayName','diffracted');
legend
title('Path Gain on a Street (AB-EW)')
xlabel('x(m)');
ylabel('Path Gain (dB)');

%plot path gain from north to south
figure; hold on;
PGCDdB = 10*log10(PGCD);
PGCD_idB = 10*log10(PGCD_i);
PGCD_ddB = 10*log10(PGCD_d);
plot(cd(2,:), PGCDdB,'DisplayName','total');
plot(cd(2,:), PGCD_idB,'DisplayName','incident');
plot(cd(2,:), PGCD_ddB,'DisplayName','diffracted');
legend
title('Path Gain on a Street (CD-NS)')
xlabel('y(m)');
ylabel('Path Gain (dB)');
