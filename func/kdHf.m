function KD=kdHf(T,par)
X=1000./T;
KD_Hf=exp(11.29e3./T  - 2.275); % true Kd_Hf from this model 2022
KD_Ti=exp(-11.05e3./T + 6.06); % Kd for Ti zircon/melt based on Ferry and Watson 2007 
KD_Y= exp(19.47 .* X - 13.04); 
KD_U= exp(15.32 .* X - 9.17); %U
KD_Th=exp(13.02e3./T -8.54); % Kd for Th
Csat=ZrSaturation(T);
KD_Sm=(13.338*Csat^(-0.622));
KD_Dy=(2460.0*Csat^(-0.867));
KD_Yb=(33460.*Csat^(-1.040));
KD_P= exp(7.646 .* X - 5.047);

switch  par.Trace
    case 'Hf'
        KD=KD_Hf;
    case 'Y'
        KD=KD_Y;
    case 'U'
        KD=KD_U;
    case 'P'
        KD=KD_P;
    case 'Sm'
        KD=KD_Sm;
    case 'Dy'
        KD=KD_Dy;
    case 'Yb'
        KD=KD_Yb;
    case 'Th'
        KD=KD_Th;
    case 'Ti'
        KD=KD_Ti;


end

end