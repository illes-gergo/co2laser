function [deff_] = deff(cry)

if cry == 0 % LN
    deff_= 168e-12;
elseif cry == 2 % ZnTe
    deff_ = 0;
elseif cry == 3 % GaP
    deff_ = 0;
elseif cry == 4 % GaAs
    deff_ = 2/sqrt(3)*42.35e-12;
elseif cry == 7 % ZnSe
    deff_ = 0;
end;

end

