function out = V5_fgv(z,A_kompozit,omega,T,k_omega,k_OMEGA,k_omegaSH,khi_eff,dnu,domega,k_omega0,omega0,gamma,cry,simp)

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
temp1 = zeros(size(ATHz));
[~,I] = max(abs(Aop));

At = ifft(Aop.*exp(-0*1i*(k_omega-k_omega0)*z*1)*2*pi*dnu*length(omega));
n2pm = fft(1i*e0*omega0*neo(2*pi*3e8/omega0,T,cry)*n2/2*abs(At).^2.*At)/dnu/2/pi/length(omega);


for nagy_omega = 2:ceil(10e12/dnu)
    temp1(nagy_omega) = -1*abszorpcio(nagy_omega)/2*ATHz(nagy_omega)-1*1i*khi_eff*omega(nagy_omega).^2/2/c^2/k_OMEGA(nagy_omega)...
        *sum(Aop(nagy_omega:end-1).*conj(Aop(1:end-nagy_omega))...
        .*exp(-1i*(k_omega(nagy_omega:end-1)-k_omega(1:end-nagy_omega)-k_OMEGA(nagy_omega))*z))*domega;
end
temp2 = zeros(size(temp1));

%for kis_omega = I-1000:I+1000
%for kis_omega = 1:end(omega)
%gamma = 0;
for kis_omega = I-simp:I+simp
%+n2pm(kis_omega)
    temp2(kis_omega) = -1*n2pm(kis_omega) -1*1i*khi_eff*omega(kis_omega).^2/2/c^2/k_omega(kis_omega)*(sum(Aop(kis_omega:end-1).*conj(ATHz(1:end-kis_omega))...
        .*exp(-1i*(k_omega(kis_omega:end-1)-k_omega(kis_omega)-k_OMEGA(1:end-kis_omega))*z))...
        +sum(Aop(kis_omega:-1:1).*(ATHz(1:kis_omega))...
        .*exp(-1i*(k_omega(kis_omega:-1:1)-k_omega(kis_omega)+k_OMEGA(1:kis_omega))*z)))*domega+...
        -1*cos(gamma)*1i*deff_*omega(kis_omega).^2/c^2/k_omega(kis_omega)*sum(ASH(kis_omega:end-1).*conj(Aop(1:end-kis_omega))...
        .*exp(-1i*(k_omegaSH(kis_omega:end-1)-k_omega(kis_omega)-k_omega(1:end-kis_omega))*z*cos(gamma)^2))*domega;
end
temp3 = zeros(size(temp1));

for kis_omega = 2*I-simp:2*I+simp
%+n2pm(kis_omega)
    temp3(kis_omega) = -1*cos(gamma)*1i*deff_*omega(kis_omega).^2/2/c^2/k_omega(kis_omega)*sum(Aop(1:kis_omega).*(Aop(kis_omega:-1:1))...
        .*exp(-1i*(k_omega(1:kis_omega)-k_omegaSH(kis_omega)+k_omega(kis_omega:-1:1))*z*cos(gamma)^2))*domega;
end

out = zeros(1,length(omega),3);
out(1,:,1) = temp1;
out(1,:,2) = temp2;
out(1,:,3) = temp3;
end