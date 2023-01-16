%LN kristály törésmutatója Sellmeier módszerrel
%0 - LN
%1 - LT
%2 - ZnTe
%3 - GaP
%4 - GaAs 
%5 - Se Edward J. Danielewicz and Paul D. Coleman, "Far Infrared Optical
%Properties of Selenium and Cadmium Telluride," Appl. Opt. 13, 1164-1170
%(1974)
%6 - CaTe


function[ er_ ] = er(omega,T,cry)

nu = omega/2/pi/3e8*0.01;

if cry == 2 %ZnTe
    
A = [4.262e-2 1.193e-2 3 0.008 0.0029 0.005];
oi = [56 94.82 182 239 299 351];
gi = [13.06 8.18 0.9 6.7 9.3 10.3];
es = 6.8;

tsc = 200e-15;
meff = 0.11*9.109e-31;
q = 1.602e-19;
e0 = 8.8541878e-12;

er_ = es+A(1)*oi(1).^2./(oi(1).^2-omega.^2*0.01^2/(2*pi*3e8)^2-2*1i*gi(1)*omega*0.01/(2*pi*3e8))+...
        A(2)*oi(2).^2./(oi(2).^2-omega.^2*0.01^2/(2*pi*3e8)^2-2*1i*gi(2)*omega*0.01/(2*pi*3e8))+...
        A(3)*oi(3).^2./(oi(3).^2-omega.^2*0.01^2/(2*pi*3e8)^2-2*1i*gi(3)*omega*0.01/(2*pi*3e8))+...
        A(4)*oi(4).^2./(oi(4).^2-omega.^2*0.01^2/(2*pi*3e8)^2-2*1i*gi(4)*omega*0.01/(2*pi*3e8))+...
        A(5)*oi(5).^2./(oi(5).^2-omega.^2*0.01^2/(2*pi*3e8)^2-2*1i*gi(5)*omega*0.01/(2*pi*3e8))+...
        A(6)*oi(6).^2./(oi(6).^2-omega.^2*0.01^2/(2*pi*3e8)^2-2*1i*gi(6)*omega*0.01/(2*pi*3e8));

er_(isnan(er_)) = 1;
er_(isinf(er_)) = 1;
end;


if cry == 4 %GaAs
    if T == 300 %ord
        e_inf = 11;
        nu_L = 292.1;
        nu_T = 268.7;
        G = 2.4;
        nu = omega/2/pi/3e8*1e-2;
        
        er_ = e_inf*(1+(nu_L^2-nu_T^2)./(nu_T^2-nu.^2+1i*G*nu));
    end;    
end;

if cry == 5
    if T == 300 %ord
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
    end;    
end;

if cry == 7
        Ein=5.9;
        rho=2.9;
        v_To=207;
        gamma0=6.50;
        v=omega/2/pi/3e8*1e-2;
        er_=Ein+(rho*v_To.^2)./(v_To.^2-v.^2+1i*gamma0*v);
end;

end

