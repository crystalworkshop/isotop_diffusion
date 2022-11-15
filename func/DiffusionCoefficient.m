function Di=DiffusionCoefficient(T,x,par)  %defining Zr diffusion coefficients in melts as f(T,X2O)
global DGfZr
theta=1000/T;
lnD=-(11.4*x+3.13)/(0.84*x+1)-(21.4*x+47)/(1.06*x+1)*theta;%best fit of Zr Diff coefficients (several workers) and WH83 dependence on XH2O
Dif=exp(lnD)*1e4*365*24*3600;% in cm2/y
mass=[89.9047026,90.9056439,91.9050386,93.9063148,95.908275];
bet=0.054+0.059;
Di(1:5)=Dif*(mass(1)./mass).^bet;
% lnHf=-3.52 - 231.09/8.31/theta;
% lnD_Hf=(-8.620340372*T-42705.17449-.318918919*x*T+4049.500765*x)/T;
Di(6)=Dif(1)*DGfZr; %exp(lnD_Hf)*1e4*365*24*3600;% in cm2/y
end