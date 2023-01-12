clc;clear;close all;
N = 1e4;

Aop = rand(1,N)+1i*rand(1,N);
k_omega = rand(1,N);
k_OMEGA = rand(1,N);
z = rand;
temp2 = conv((Aop),(Aop),"full");
temp1 = zeros(1,N);

% for kis_omega =1:length(Aop)
%     temp1(kis_omega) = sum(conj(Aop(kis_omega:-1:1)).*(Aop(1:kis_omega)));
% end
% temp2 = (conv(conj((Aop)),((Aop)),"full"));
for nagy_omega = 1:length(Aop)
    temp1(nagy_omega) = sum(Aop(nagy_omega:end-1).*conj(Aop(1:end-nagy_omega))...
         .*exp(-1i*(k_omega(nagy_omega:end-1)-k_omega(1:end-nagy_omega)-k_OMEGA(nagy_omega))*z));
                            
end
temp2 = flip(conv(conj(Aop).*exp(-1i.*k_omega.*z),(flip(Aop.*exp(1i*k_omega.*z))),"full"));
temp2 = [temp2(N:end)].*exp(1i.*k_OMEGA.*z);

plot(real(temp1))
hold on
plot(real(temp2))
plot(real(temp2-temp1))
hold off;
figure
plot(imag(temp1))
hold on
plot(imag(temp2))
plot(imag(temp2-temp1))