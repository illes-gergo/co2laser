clc;clear;close all;delete(gcp('nocreate'));
%parpool("threads");

%z_ARRAY = [2,4,8]*1e-3;
%for ITER = 1:length(z_ARRAY)
%    clearvars("-except","ITER","z_ARRAY")
tic;
cases = [[2e-3;100e13],[4e-3;100e13],[4e-3;50e13],[4e-3;25e13],[8e-3;25e13]];
parfor ITER = 1:length(cases)
%adatok kiírása file-ba
dir_n = convertCharsToStrings(strrep(strrep(datestr(datetime), ' ', '_'),':','.'));
dir_n = num2str(ITER)+dir_n;
mkdir(dir_n);
cry = 4; % 4 - GaAs  7 - ZnSe  2 - ZnTe

T = 300;    %K
c = 3e8;    %m/s
c0 = 3e8;
lambda0 = 10.6e-6;   %m
N = 4*1e4;    %db
tau = 1.75e-12;  %s
I0 = cases (2,ITER);%100e13;%100/sqrt(tau/100e-15)*1e13;    %GW/cm^2
%I0 = 100e13;
khi_eff =   2*deffTHz(cry);%78.4e-12;%2*65.6e-12;%360e-12; %pm/V;
e0 = 8.854e-12;  %F*m^2
nu0 = 1.5e12;
deltanu = nu0; %1e12;
simp = 1100;


dz = 0.5e-5;
%z_vegso = 4e-3;
z_vegso= cases(1,ITER);%8e-3;
z = 0:dz:z_vegso;
omega0 = 2*pi*c/lambda0;

elochirp = 1*z_vegso/2;

utem = fix(length(z)/10);
if utem == 0 
    utem = 1;
end


omegaMAX = 5*2*omega0;%5e14*2*pi;
dt = 2*pi/omegaMAX;
t = (0:N-1)*dt;
t = t-t(end)/2;
domega = omegaMAX/N;
dnu = domega/2/pi;

omega = (0:N-1)*domega;
nu=omega/2/pi;

deltaOmega =2*sqrt(2*log(2))/tau;

lambda = 2*pi*c./omega;
lambda(1) = lambda(2);
ngp0 = ngp(lambda0,T,cry);
np0 = neo(lambda0,T,cry);
ngpSH = ngp(lambda0/2,T,cry);
npSH = neo(lambda0/2,T,cry);

nTHz = nTHzo(2*pi*nu0,T,cry);
vfTHz = c./nTHz;

gamma = acos(ngp0/nTHzo(2*pi*nu0,T,cry));

fileID = fopen(strcat(dir_n,'/data.txt'),'w');
fprintf(fileID,'CO2 lézerben történő THz-keltés GaAs-ben\n\n');
fprintf(fileID,'Pumpa hullámhossza: %6.2f um\n',lambda0*1e6);
fprintf(fileID,'Pumpa intenzitása: %6.2f GW/cm^2\n',I0*1e-13);
fprintf(fileID,'Pumpa FL félértékszélessége: %6.2f ps\n',tau*1e12);
fprintf(fileID,'Fázisillesztési frekvencia (impulzusfrontdöntéssel): %6.2f THz\n',nu0*1e-12);
fprintf(fileID,'Impulzusfrontdöntés: %6.2f °\n\n',gamma/pi*180);

fprintf(fileID,'Kristály anyaga: GaAs\n');
fprintf(fileID,'Kristályhossz: %6.2f mm\n',z_vegso*1e3);
fprintf(fileID,'Kristály hőmérséklete: %3.0f\n',T);
fprintf(fileID,'Pumpa fázis törésmutatója %6.2f nm-es hullámhosszon: %6.4f\n',lambda0*1e9,np0);
fprintf(fileID,'Pumpa csoport törésmutatója %6.2f nm-es hullámhosszon: %6.4f\n',lambda0*1e9,ngp0);
fprintf(fileID,'Másodharmonikus fázis törésmutatója %6.2f nm-es hullámhosszon: %6.4f\n',lambda0/2*1e9,npSH);
fprintf(fileID,'Másodharmonikus csoport törésmutatója %6.2f nm-es hullámhosszon: %6.4f\n',lambda0/2*1e9,ngpSH);

fprintf(fileID,'THz fázis törésmutatója %6.2f THz-es frekvencián: %6.4f\n',nu0*1e-12,nTHz);
fprintf(fileID,'Nemlineáris szuszceptibilitás (Chi_eff): %6.2f pm/V\n\n',khi_eff*1e12);
fprintf(fileID,'Nemlineáris együttható (d_eff): %6.2f pm/V\n\n',deff(cry)*1e12);
fprintf(fileID,'Nemlineáris törésmutató (n2): %6.8f cm^2/GW\n\n',n2value(cry)*1e13);

fprintf(fileID,'Térbeli lépések száma: %d (dz = %f6.2 um)\n',length(z), dz*1e6);
fprintf(fileID,'Spektrális felbontás: %d (d_nu = %f6.2 GHz)\n',N, dnu*1e-9);
fprintf(fileID,'Spektrális tartomány: (Nu_max=) %6.2 THz',omegaMAX/2/pi*1e-12);
fprintf(fileID,'Időbeli felbontás: %d (d_t = %f6.2 fs)\n',N, dt*1e15);
fprintf(fileID,'Időbeli tartomány: (T_max=) %6.2 ps\n\n',N*dt*1e12);

fprintf(fileID,'Számítás időpontja: %s\n',dir_n);

%Adatok kiírása file-ba

%A0 = sqrt(2*I0/neo(lambda0,T,cry)/e0/c)*tau/(2*sqrt(2*pi*log(2)));
A0t = sqrt(2*I0/neo(lambda0,T,cry)/e0/c);
%Aot = A0t*exp(-2*log(2)*t.^2/tau^2).*exp(1i*(omega0)*t);
%Aop = fft(Aot)*dt/2/pi;
%Aop = Aop.*exp(1i*2*2200e-30/2*(omega-omega0).^2);
%Aop = A0*exp(-((omega-omega0).^2/deltaOmega.^2));
%Aop0 = Aop;

%gamma = 0;
n_omega = neo(lambda,T,cry);
k_OMEGA = real(omega.*nTHzo(omega,T,cry)/c);%+1e5;
k_OMEGA0 = real(omega.*nTHzo(2*pi*nu0,T,cry)/c);
ddk_omega = -ngp0.^2/omega0/c/np0*tan(gamma)^2;
k_omega = real(1/cos(gamma).*(omega.*n_omega/c+1*(omega-omega0).^2/2.*ddk_omega));%+1e5;
ddk_omegaSH = -ngpSH.^2/omega0/2/c/npSH*tan(gamma)^2;
k_omegaSH = real(1/cos(gamma).*(omega.*n_omega/c+1*(omega-2*omega0).^2/2.*ddk_omegaSH));%+1e5;
k_omega0 = real(1/cos(gamma).*(omega.*ngp0/c));
k_omegaSH0 = real(1/cos(gamma).*(omega.*ngpSH/c));

ko0 = interp1(omega,k_omega,omega0);
ksh = interp1(omega,k_omega,2*omega0);
dk2 = 2*k_omega(1:end/2)-k_omegaSH(1:2:end);
dk3 = 2*neo(lambda,T,cry).*omega/c-2*neo(lambda/2,T,cry).*omega/c;
%dk = 2*real(1/cos(gamma).*(omega.*neo(lambda,T,cry)/c+(omega-omega0).^2/2.*ddk_omega))...
%    -real(1/cos(gamma).*(2*omega.*neo(lambda/2,T,cry)/c+(2*omega-2*omega0).^2/2.*ddk_omegaSH));
%return;
%k_omegaSH = real(1/cos(gamma).*(omega.*ngp0/c));


A0 = sqrt(2*I0/neo(lambda0,T,cry)/e0/c)*tau/(2*sqrt(2*pi*log(2)));
%Aop = A0*exp(-((omega-omega0).^2/deltaOmega.^2)).*exp(1i*(k_omega-k_omega0)*elochirp);
%Aop0 = Aop;
Aot = A0t*exp(-2*log(2)*t.^2/tau^2).*exp(1i*omega0*t);
Aop = fft(Aot).*exp(1i*(k_omega-k_omega0)*elochirp)*dt/2/pi;
Aop0 = Aop;

%for j = 1:10
%for kis_omega = 7000:9000
%+n2pm(kis_omega)
%    temp3(kis_omega,j) = sum(Aop(1:kis_omega).*(Aop(kis_omega:-1:1)).*...
%        exp(-1i*(k_omegaSH(1:kis_omega)-k_omegaSH(kis_omega)+k_omega(kis_omega:-1:1))*z(50)))*domega;
%end;
%end;
%return;


ATHz = zeros(size(Aop));
ASH = zeros(size(Aop));

pF = sum(abs(Aop).^2)*np0;
%FISH = neo(lambda,T,cry);
FI = 4*nTHz.^2./(1+nTHz).^2;
FA = 2*nTHz./(1+nTHz);
%return;

    dlmwrite(strcat(dir_n,'/PumpInt0.txt'),[t.' np0*e0*c0/2*abs((ifft(Aop)).').^2*omegaMAX^2]);
    dlmwrite(strcat(dir_n,'/PumpSpec0.txt'),[lambda.' (abs(Aop).').^2]);
    
A_komp(1,:,1) = ATHz;
A_komp(1,:,2) = Aop0;
A_komp(1,:,3) = ASH;

%% kezdĹ komplex tĂŠrerĹssĂŠgek ĂśsszeillesztĂŠse
%% diffegyenlet lĂŠtrehozĂĄsa
v6_fgv =@(z,A_kompozit) diffegy_conv(z,A_kompozit,omega,T,k_omega,k_OMEGA,k_omegaSH,khi_eff,dnu,domega,k_omega0,omega0,gamma,cry,simp);

%% diffegyenlet megoldĂĄsa

subplot(2,4,3);

%plot(1e-12*nu,nTHzo(2*pi*nu,T,cry));
%xlim([0 5]);
%title('THz refractive index');

subplot(2,4,4);

%plot(1e-12*nu,1e-2*aTHzo(2*pi*nu,T,cry));
%xlim([0 5]);
%title('THz absorption (1/cm)');
effic = zeros(size(z));
efficSH = zeros(size(z));
for ii = 1:length(z)
A_komp(1,:,1) = ATHz;
A_komp(1,:,2) = Aop;
A_komp(1,:,3) = ASH;
% tic;
[z2, A_komp] = RK4_M(v6_fgv,dz,(ii-1)*dz,A_komp,(ii+0.1)*dz);
% toc;
ATHz = A_komp(2,:,1).';
Aop = A_komp(2,:,2).';
ASH = A_komp(2,:,3).';
%return;
A_komp = zeros(1,length(omega),3);

effic(ii) = sum(abs(ATHz).^2.*FI.')/pF ;
efficSH(ii) = sum(abs(ASH).^2*npSH.')/pF ;
%% kirajzolĂĄs

subplot(2,4,1)
plot(lambda,abs(Aop).^2);
xlim([(lambda0-1500e-9),lambda0+1500e-9]);
title('Pump Intensity spectrum');
xlabel('Hullámhossz (m)');

subplot(2,4,2)
plot(t*1e12, 1e-13*np0*e0*c0/2*abs((ifft(Aop.*exp(-1i*(k_omega-k_omega0)*z(ii)).'))*omegaMAX).^2);
xlim([-5 5]);
%xlim([omega0-omega0/10,omega0+omega0/10]/2/pi);
title('Pump Intensity');
xlabel('Time (ps)');
ylabel('Intensity (GW/cm^2)');
%return;

subplot(2,4,7);
plot(1e-12*nu,abs(ATHz));
xlim([0,5]);
title('THz amplitude spectrum');
xlabel('Frequency (THz)');

subplot(2,4,5);
b = plot(z(1:ii)*1e3,effic(1:ii)*100);
title('Efficiency (%)');
xlabel('Crystal length (mm)');


subplot(2,4,6);
b = plot(z(1:ii)*1e3,efficSH(1:ii)*100);
title('SH-Efficiency (%)');
xlabel('Crystal length (mm)');

subplot(2,4,8)
plot(t*1e12,1e-5*real((ifft(FA.'.*ATHz.*exp(-1i*(k_OMEGA-k_OMEGA0)*z(ii)).')))*omegaMAX);
title('THz electric field (inside)');
xlabel('Time (ps)');
ylabel('Electric Field (KV/cm)');
xlim([-10 10]);


subplot(2,4,3);
plot(lambda,abs(ASH).^2);
xlim([(lambda0/2-1000e-9),lambda0/2+1000e-9]);
title('SH Intensity spectrum');
xlabel('Hullámhossz (m)');

%plot(1e-12*nu,nTHzo(2*pi*nu,T,cry));
%xlim([0 5]);
%title('THz refractive index');

subplot(2,4,4);
plot(t*1e12, 1e-13*neo(lambda0/2,T,cry)*e0*c0/2*abs((ifft(ASH.*exp(-1i*(k_omegaSH-k_omegaSH0)*z(ii)).'))*omegaMAX).^2);
xlim([-5 5]);
%xlim([omega0-omega0/10,omega0+omega0/10]/2/pi);
title('SH Intensity');
xlabel('Time (ps)');
ylabel('Intensity (GW/cm^2)');


drawnow;

if mod(ii,utem)==0
    dlmwrite(strcat(dir_n,'/PumpInt-',int2str(fix(ii*dz*1e6)),'_um.txt'),[(t).' np0*e0*c0/2*abs((ifft(Aop.*exp(-1i*(k_omega-k_omega0)*z(ii)).'))*omegaMAX).^2]);
    dlmwrite(strcat(dir_n,'/PumpSpec-',int2str(fix(ii*dz*1e6)),'_um.txt'),[lambda.' abs(Aop).^2]);
    dlmwrite(strcat(dir_n,'/THzt-',int2str(fix(ii*dz*1e6)),'_um.txt'),[t.' real((ifft(FA.'.*ATHz.*exp(-1i*(k_OMEGA-k_OMEGA0)*z(ii)).')))*omegaMAX]);
    dlmwrite(strcat(dir_n,'/THzSpec-',int2str(fix(ii*dz*1e6)),'_um.txt'),[nu.' abs(ATHz)]);
end;

end

hossz = toc;


h = fix(hossz/3600);
min = fix((hossz-h*3600)/60);
sec = hossz-h*3600-min*60;
fileID = fopen(strcat(dir_n,'/data.txt'),'a+');
fprintf(fileID,'/nSzámítás hossza: %2.0f óra %2.0f perc %2.0f másodperc', h, min, sec);
fclose(fileID);


    dlmwrite(strcat(dir_n,'/PumpInt.txt'),[(t).' np0*e0*c0/2*abs((ifft(Aop.*exp(-1i*(k_omega-k_omega0)*z(ii)).'))).^2*omegaMAX^2]);
    dlmwrite(strcat(dir_n,'/PumpSpec.txt'),[lambda.' abs(Aop).^2]);
    dlmwrite(strcat(dir_n,'/THzt.txt'),[t.' real((ifft(FA.'.*ATHz.*exp(-1i*(k_OMEGA-k_OMEGA0)*z(ii)).')))*omegaMAX]);
    dlmwrite(strcat(dir_n,'/THzSpec.txt'),[nu.' abs(ATHz)]);
    dlmwrite(strcat(dir_n,'/efficTHz.txt'),[z.' effic.']);
    dlmwrite(strcat(dir_n,'/efficSH.txt'),[z.' efficSH.']);


end
