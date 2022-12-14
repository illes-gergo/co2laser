
%%
clc;clear;close all;
N = 1E4;
Aop1 = rand(1,N)+1i*rand(1,N);
Aop2 = rand(1,N)+1i*rand(1,N);
eredmeny1 = zeros(1,N);
for ii = 1:length(Aop1)
    eredmeny1(ii) = sum(Aop1(ii:-1:1).*Aop1(1:ii));
end
eredmeny2 = conv(Aop1,Aop1);
nexttile;
plot(abs(eredmeny1))
nexttile;
plot(abs(eredmeny2));
nexttile;
plot(imag(eredmeny1-eredmeny2(1:N))./imag(eredmeny1))
%%
clc;clear;close all;
N = 1E4;
Aop1 = rand(1,N)+1i*rand(1,N);
Aop2 = rand(1,N)+1i*rand(1,N);
eredmeny1 = zeros(1,N);
for ii = 1:length(Aop1)
    eredmeny1(ii) = sum(Aop1(ii:end).*conj(Aop2(1:end-ii+1)));
end
eredmeny2 = conv(conj(flip(Aop2)),(Aop1));
nexttile;
plot(abs(eredmeny1))
nexttile;
plot(abs(eredmeny2));
nexttile;
plot(imag(eredmeny1-eredmeny2(N:end))./imag(eredmeny1))
nexttile;
plot(abs(eredmeny1))
hold on;
plot(abs(eredmeny2(N:end)))

%%
clc;clear;close all;delete(gcp('nocreate'));
N = 1E4;
Aop1 = rand(1,N)+1i*rand(1,N);
Aop2 = rand(1,N)+1i*rand(1,N);
k_omega = rand(1,N);
k_OMEGA = rand(1,N);
eredmeny1 = zeros(1,N);
z = rand*100;
tic;
for ii = 1:length(Aop1)
    eredmeny1(ii) = sum(Aop1(ii:end)...
        .*conj(Aop1(1:end-ii+1))...
        .*exp(-1i*(k_omega(ii:end)-k_omega(1:end-ii+1)-k_OMEGA(ii))*z));
end
toc;
tic;
eredmeny2 = conv(flip(conj(Aop1).*exp(1i.*k_omega.*z)),(Aop1.*exp(-1i*k_omega.*z)),"full");
eredmeny2 = eredmeny2(N:end).*exp(1i.*k_OMEGA.*z);
toc;
nexttile;
plot(abs(eredmeny1));
nexttile;
plot(real(eredmeny1))
hold on
plot(real(eredmeny2));
nexttile;
plot(real(eredmeny1-eredmeny2)./real(eredmeny1))
nexttile;
plot(imag(eredmeny1-eredmeny2)./imag(eredmeny1))

