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

function [ nTHz ] = nTHzo( omega, T,cry)

if cry == 0

    if T == 100    
        nTHz = 4.69732-0.03006e-12*omega/2/pi+0.04066e-24*(omega/2/pi).^2;
    elseif T == 300
        nTHz = 4.91372-0.01747e-12*omega/2/pi+0.04004e-24*(omega/2/pi).^2;
    end;

elseif cry == 2
    
    nTHz = real(sqrt(er(omega,T,cry)));
    
elseif cry == 4
    
    nTHz = real(sqrt(er(omega,T,cry)));
    
elseif cry == 7
    
    nTHz = real(sqrt(er(omega,T,cry)));
    
end;
    
end

