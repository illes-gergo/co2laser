function [n2_] = n2value(cry)

if cry == 0 % LN
    n2_= 5.3e-19;
elseif cry == 2 % ZnTe
    n2_ = 9e-18; % approx
%elseif cry == 3 % GaP
    %deff_ = 24.8e-12;
elseif cry == 4 % GaAs
    n2_ = 5.9e-18; %10.6 um-en  %1.7e-17;%42.35e-12;%65.6e-12
elseif cry == 7 % ZnSe
    n2_ = 4.7e-18;%73e-12;
end;

end