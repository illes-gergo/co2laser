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

function [ alpha ] = aTHzo( omega, T,cry)
if cry == 0
    if T == 100
        [a b] = min(abs(omega-1.5708e13));
        alpha(1:b) = 100*(6.20304e-12*omega(1:b)/2/pi-5.03352e-24*(omega(1:b)/2/pi).^2+2.63295e-36*(omega(1:b)/2/pi).^3);
        alpha(b:length(omega)) = 100*(6.20304e-12*omega(b)/2/pi-5.03352e-24*(omega(b)/2/pi).^2+2.63295e-36*(omega(b)/2/pi).^3);    
    elseif T == 300
        alpha = 100*(2.16411e-12*omega/2/pi+10.81e-24*(omega/2/pi).^2);
    end;
    
elseif cry == 2
    alpha = 2*omega/3e8.*imag(sqrt(er(omega,T,cry)));

    
elseif cry == 4
    alpha = -2*omega/3e8.*imag(sqrt(er(omega,T,cry)));

elseif cry == 7
    alpha=-2*omega/3e8.*imag(sqrt(er(omega,T,cry)));
end;

end

