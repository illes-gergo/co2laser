function out = diffegy_conv(z,A_kompozit,omega,T,k_omega,k_OMEGA,k_omegaSH,khi_eff,dnu,domega,k_omega0,omega0,gamma,cry,simp)

% T = 100;    %K
n2 = n2value(cry);
deff_ = 1*deff(cry);
c = 3e8;    %m/s
e0 = 8.854e-12;
% lambda0 = 1031.8e-9;   %m
% N = 2*1e4;    %db
% tau = 150e-15;  %s
% I0 = 20e1/tau;    %GW/cm^2
% khi_eff =  360e-12; %pm/V;
% e0 = 8.854e-12;  %F*m^2
% nu0 = 0.5e12;
%
% omega0 = 2*pi*c/lambda0;
% omegaMAX = 5e14*2*pi;
% domega = omegaMAX/N;
% dnu = domega/2/pi;
%
% omega = (0:N-1)*domega;
%
% deltaOmega =2*sqrt(2*log(2))/tau;
%
% lambda = 2*pi*c./omega;
% lambda(1) = lambda(2);
% ngp0 = ngp(lambda0,T);
% np0 = neo(lambda0,T);
%
% nTHz = nTHzo(lambda0,T);
% vfTHz = c./nTHz;
%
% gamma = acos(ngp0/nTHzo(2*pi*nu0,T));
%
% A0 = sqrt(2*I0/neo(lambda0,T)/e0/c)*tau*sqrt(pi/log(2)); % ???
% A0 = sqrt(2*I0/neo(lambda0,T)/e0/c)*tau/(2*sqrt(2*pi*log(2)));
% Aop = A0*exp(-((omega-omega0).^2/deltaOmega.^2));
% plot(Aop(:,1));
% return;
% n_omega = neo(lambda,T);
% k_OMEGA = real(omega.*nTHzo(omega,T)/c);%+1e5;
% ddk_omega = -ngp0.^2/omega0/c/np0*tan(gamma)^2;
% k_omega = real(1/cos(gamma).*(omega.*n_omega/c+(omega-omega0).^2/2.*ddk_omega));%+1e5;
%
%
%
% if z_index == 1
%

% end

abszorpcio = aTHzo(omega,T,cry);
abszorpcio(abszorpcio>1e5) = 1e5;

ATHz = A_kompozit(1,:,1);
Aop = A_kompozit(1,:,2);
ASH = A_kompozit(1,:,3);
%temp1 = zeros(size(ATHz));
[~,I] = max(abs(Aop));
NN = length(Aop);
At = ifft(Aop.*exp(-0*1i*(k_omega-k_omega0)*z*1)*2*pi*dnu*length(omega));
n2pm = fft(1i*e0*omega0*neo(2*pi*3e8/omega0,T,cry)*n2/2*abs(At).^2.*At)/dnu/2/pi/length(omega);
% 
% for nagy_omega = 2:ceil(10e12/dnu)
%     temp1(nagy_omega) = -1*abszorpcio(nagy_omega)/2*ATHz(nagy_omega)-1*1i*khi_eff*omega(nagy_omega).^2/2/c^2/k_OMEGA(nagy_omega)...
%         .*sum(Aop(nagy_omega:end).*conj(Aop(1:end-nagy_omega+1))...
%         .*exp(-1i*(k_omega(nagy_omega:end)-k_omega(1:end-nagy_omega+1)-k_OMEGA(nagy_omega))*z))*domega;
% end
%tic;
temp11 = conv(flip(conj(Aop).*exp(1i.*k_omega.*z)),(Aop.*exp(-1i*k_omega.*z)),"full");
temp11 = temp11(NN:end).*exp(1i.*k_OMEGA.*z).*(-1.*1i.*khi_eff.*omega.^2/2/c^2./k_OMEGA).*domega-1.*abszorpcio/2.*ATHz;
temp11(1) = 0;
temp1 = temp11;
%toc;
% close all;
% CHECK = temp1-temp11;
% subplot(2,2,1);
% plot(real(CHECK)./real(temp1));
% subplot(2,2,2);
% plot(imag(CHECK)./imag(temp1));
% subplot(2,2,3)
% plot(real(temp1));
% hold on
% plot(real(temp11));
% hold off;
% subplot(2,2,4)
% plot(imag(temp1));
% hold on
% plot(imag(temp11));
% hold off;
% 
% drawnow;
% tic;
temp2 = zeros(size(temp1));
%
% %for kis_omega = I-1000:I+1000
% %for kis_omega = 1:end(omega)
% %gamma = 0;
% for kis_omega = I-simp:I+simp
% %+n2pm(kis_omega)
%     temp2(kis_omega) = -1*n2pm(kis_omega) -1*1i*khi_eff*omega(kis_omega).^2/2/c^2/k_omega(kis_omega)*(sum(Aop(kis_omega:end-1).*conj(ATHz(1:end-kis_omega))...
%         .*exp(-1i*(k_omega(kis_omega:end-1)-k_omega(kis_omega)-k_OMEGA(1:end-kis_omega))*z))...
%         +sum(Aop(kis_omega:-1:1).*(ATHz(1:kis_omega))...
%         .*exp(-1i*(k_omega(kis_omega:-1:1)-k_omega(kis_omega)+k_OMEGA(1:kis_omega))*z)))*domega+...
%         -1*cos(gamma)*1i*deff_*omega(kis_omega).^2/c^2/k_omega(kis_omega)*sum(ASH(kis_omega:end-1).*conj(Aop(1:end-kis_omega))...
%         .*exp(-1i*(k_omegaSH(kis_omega:end-1)-k_omega(kis_omega)-k_omega(1:end-kis_omega))*z*cos(gamma)^2))*domega;
% end
% toc;
%tic;
temp21 = conv(flip(conj(ATHz).*exp(-1i.*k_OMEGA.*z)),Aop.*exp(-1i.*k_omega.*z));
temp21 = temp21(NN:end).*exp(1i.*k_omega.*z);
temp22 = conv(Aop.*exp(-1i.*k_omega.*z),ATHz.*exp(-1i.*k_OMEGA.*z));
temp22 = temp22(1:NN).*exp(1i.*k_omega.*z);
temp20 = -1*n2pm -1*1i*khi_eff.*omega.^2/2/c^2./k_omega.*(temp21+temp22).*domega;
%temp20(1) = 0;
temp23 = conv(flip(conj(Aop).*exp(1i.*k_omega.*z.*cos(gamma).^2)),ASH.*exp(-1i.*k_omegaSH.*z.*cos(gamma).^2)).*domega;
temp23 = -1*cos(gamma).*1i.*deff_.*omega.^2/c^2./k_omega.*temp23(NN:end).*exp(1i.*k_omega.*z.*cos(gamma).^2);
temp24 = temp20+temp23;
temp24(1) = 0;
temp2 = temp24;
% close all;
% CHECK = temp2-temp24;
% subplot(2,2,1);
% plot(real(CHECK)./real(temp2));
% subplot(2,2,2);
% plot(imag(CHECK)./imag(temp2));
% subplot(2,2,3)
% plot(real(temp2));
% hold on
% plot(real(temp24));
% hold off;
% subplot(2,2,4)
% plot(imag(temp2));
% hold on
% plot(imag(temp24));
% hold off;
% 
% drawnow;
% 
%toc;
% temp3 = zeros(size(temp1));
% 
% for kis_omega = 2*I-simp:2*I+simp
% %+n2pm(kis_omega)
%     temp3(kis_omega) = -1*cos(gamma)*1i*deff_*omega(kis_omega).^2/2/c^2/k_omega(kis_omega)*sum(Aop(1:kis_omega).*(Aop(kis_omega:-1:1))...
%         .*exp(-1i*(k_omega(1:kis_omega)-k_omegaSH(kis_omega)+k_omega(kis_omega:-1:1))*z*cos(gamma)^2))*domega;
% end
% tic;
temp31 = conv(Aop.*exp(-1i.*k_omega.*z*cos(gamma)^2),Aop.*exp(-1i.*k_omega.*z*cos(gamma)^2))*domega;
temp31 =  -1*cos(gamma)*1i*deff_.*omega.^2/2/c^2./k_omega.*temp31(1:NN).*exp(1i.*k_omegaSH.*z.*cos(gamma).^2);
temp31(1) = 0;
temp3 = temp31;
%toc;

% close all;
% CHECK = temp3-temp31;
% subplot(2,2,1);
% plot(real(CHECK)./real(temp3));
% subplot(2,2,2);
% plot(imag(CHECK)./imag(temp3));
% subplot(2,2,3)
% plot(real(temp3));
% hold on
% plot(real(temp31));
% hold off;
% subplot(2,2,4)
% plot(imag(temp3));
% hold on
% plot(imag(temp31));
% hold off;
% 
% drawnow;
% 

out = zeros(1,length(omega),3);
out(1,:,1) = temp1;
out(1,:,2) = temp2;
out(1,:,3) = temp3;

end

% temp1 = flip(conv(conj(Aop).*exp(1i.*k_omega.*z),(flip(Aop.*exp(1i*k_omega.*z))),"full"));
% temp1 = [temp1(NN:end)].*exp(1i.*k_OMEGA.*z).*(-1.*abszorpcio/2.*ATHz-1.*1i.*khi_eff.*omega.^2/2/c^2./k_OMEGA).*domega;