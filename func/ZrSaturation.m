function Csat=ZrSaturation(T) %defining Zr saturation conditions
%Csat=4.414e7/exp(13352/T)/2; %Watson 96, Eq 1, in ppm Zr for checking. (divide by 1),or mol Zr (divide by 2)
%Mfactor = 0.0000048*(T)^2-0.0083626*(T)+4.8484463; % empirical relations from magma Fig.
%differentiation calc (file  M_factorsforOleg.xlsx
Mfactor=1.62;
Csat=490000/exp(10108/T+1.16*(Mfactor-1)-1.48); % Boehnkeetal2013ChemGeol351,324 assuming 490000ppm in Zircon

% Mfactor=1.3;
% Csat=490000/exp(10108/T+1.16*(Mfactor-1)-1.48); % Boehnkeetal2013ChemGeol351,324 assuming 490,000 ppm in Zircon
%Csat=490000/(exp(12900/T-0.85*(Mfactor-1)-3.80));% Watson and Harrison 1983
% for Monazite (Does not work for some reason):
% H2O = 1wt% use below expression (Table 4, Rapp Watson 86)
%Csat=0.0000190*exp(0.0143872*T);
%Csat=600000/(exp(-0.0144*T+24.177));
% H2O = 6wt% use below expression (Table 4, Rapp Watson 86)
%Csat=0.00012*exp(0.01352*T);
%Csat=600000/(exp(-0.0135*T+22.296));


%for Apatite:
%SiO2=0.68;
%Csat=430000 / exp((-4800 + 26400 * SiO2) / T + 3.10 - 12.4 * SiO2);%Harrison Watson 84

end
