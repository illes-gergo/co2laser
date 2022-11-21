%LN kristály törésmutatója Sellmeier módszerrel
%0 - LN
%1 - LT
%2 - ZnTe
%3 - GaP
%4 - GaAs T. Skauli, P. Kuo, K.L. Vodopyanov eta l., J. Appl. Phys. 94,
%6447-6455 (2003)
%5 - Se Edward J. Danielewicz and Paul D. Coleman, "Far Infrared Optical
%Properties of Selenium and Cadmium Telluride," Appl. Opt. 13, 1164-1170
%(1974)

function [ ng ] = ngp( lambda,T, cry )
lambda1 = lambda*1e6;
if cry == 0
syms x;
%lambda1 = lambda*1e6;
a1 = 5.078;
a2 = 0.0964;
a3 = 0.2065;
a4 = 61.16;
a5 = 10.55;
a6 = 1.59e-2;
b1 = 4.677e-7;
b2 = 7.822e-8;
b3 = -2.653e-8;
b4 = 1.096e-4;
T = T-273;
f = (T-24.5)*(T+570.82);

n0 = sqrt(a1+b1*f+(a2+b2*f)./(x.^2-(a3+b3*f)^2)+(a4+b4*f)./(x.^2-a5^2)-a6*x.^2);
        
        
        a = n0-x*diff(n0);
        x = lambda1;
        ng = eval(a);
        
elseif cry == 2
syms l;
lambda1 = lambda*1e6;
%l = lambda*1e6;
a0 = 9.92;
a1 = 0.4253;
b1 = 0.37766;
a2 = 2.6358;
b2 = 56.5;
%H. H. Li, doi.org/10.1063/1.555705;

n0 = real(sqrt(a0+a1./(l.^2-b1^1)+a2./(l.^2/b2^2-1)));
%lambda= linspace(0.8,11,N); %wavelength in um
%n0=sqrt(A+B*l^2./(l.^2-C.^2)+D*l.^2./(l.^2-E.^2));% optical refractive index
        
        
        a = n0-l*diff(n0);
        l = lambda1;
        ng = eval(a);
        
elseif cry == 4
    syms l;
    
    a0 = 4.372514;
    a = [5.466742 0.0242996 1.957522];
    b = [0.4431307 0.8746453 36.9166];
    
    n0 = real(sqrt(a0+1+a(1)*l.^2./(l.^2-b(1)^2)+a(2)*l.^2./(l.^2-b(2)^2)+...
        a(3)*l.^2./(l.^2-b(3)^2)));
    
        a = n0-l*diff(n0);
        l = lambda1;
        ng = eval(a);
    
elseif cry == 5
        
        syms nu;
        
        nu0 = 1./lambda*0.01;
        
        Se1 = 5.79;
        Se2 = 0.18;
        So1 = 0.01;
        So2 = 0.45;
        k_inf = 6.97;
        kDC = 7.43;
        nue1 = 30648.7;
        nue2 = 12904.8;
        nuo1 = 231;
        nuo2 = 140;
        g01 = 9;
        g02 = 8;
        
        er_ = 1+Se1*nue1^2./(nue1^2-nu.^2)+Se2*nue2^2./(nue2^2-nu.^2)+...
            So1*nuo1^2./(nuo1^2-nu.^2+1i*g01*nu)+So2*nuo2^2./(nuo2^2-nu.^2+1i*g02*nu);
        
        n0 = sqrt(real(er_));
        
        a = n0+nu*diff(n0);
        nu = nu0;
        ng = eval(a);
        
elseif cry == 7
    syms l;
    A=1;
    B=4.2980149;
    C=0.1920630;
    D=0.62776557;
    E=0.3787826;
    F=2.8955633;
    G=46.994595;    
    
    n0=real(sqrt(A+B*l.^2./(l.^2-C.^2)+D*l.^2./(l.^2-E.^2)+F*l.^2./(l.^2-G.^2)));
    
        a = n0-l*diff(n0);
        l = lambda1;
        ng = eval(a);    
    end;
end

