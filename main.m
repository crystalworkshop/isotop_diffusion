clear all
close all
oldpath = path;
currentFolder = pwd;
funcpas=genpath([currentFolder,'/func']);
addpath(funcpas)
tmpName = tempname('/Results');
global  n R A B D F alpha1 beta Xs UCR ZrPl Tsat  CbulkZr MinCore DGfZr S0
global Dplag Temp MeltFrac time XH2O Tsolidus Csupsat V Dscale tscale length S W t CZirc CPl ZircNuc Czl Czh
%%parameters for simulations
CbulkZr = 100; %bulk rock ppm, set up by user
tyear=3600*24*365;
iplot=1; %plot results
n=500; %number of points in radial mesh. Can be chaged by user depends on desired accuracy
nt=500;
CZirc=490000.; %zirconium concentration in zircon, ppm
XH2O=2; %initial water content in melt, required for diffusion coefficient simulations.
Tsolidus=400+273; %arbitrary solidus temperature for phase diagram used
Csupsat=3; %ppm supersaturation to cause nucleation of a new crystal upon cooling
% Dplag=0.4; %Distribution coefficient for main minerals (0 -fully incompatoble)
UCR=1; %Critical concentration for microZircon nucleation on major minerals
ZircNuc=1e-4; % Ziron stable nuclei in cm,
length=0.1; %20e-4*(CZirc/CbulkZr)^(1./3.); radius of melt cell
DGfZr=0.5; %ratio of diffusion coefficients of Hf to Zr; change for other element of interest

fig_th=figure('Units','normalized','Color','w','Position' ,[0.05,0.05,0.3,0.85]);
%options = optimoptions('fsolve','Display','off');
Tsat=fsolve(@(T)ZrSaturation(T)*mf_rock(T-273.15)-CbulkZr,1000);

% parameters for the simuklation
par.Tsat=Tsat; %Starting at saturation
par.Tend=695;  %final temerature,C
par.tfin=1500; %final time
par.Cbulk=CbulkZr;
par.RhoZrM=4.7/2.3; %Ratio of zircon to melt density
par.Kmin=0.1;  %Parameters for zirconiun partition coefficient in magor phase
par.Kmax=0.1;
par.Crit=30;
par.delta=0.2;
par.Ktrace=0.1; %trace partition coefficient in magor phase.
par.Trace='Hf';
par.XH20=XH2O;
par.L=length;
par.DGfZr=DGfZr;  %diffusion coefficient ratio
par.nt=nt;

tr=[];
Tr=[];
[time, Temp,MeltFrac]=TemperatureHistory_m_erupt(tr,Tr,par);

%% :Scaling
tfin=time(end); %total time in years of the process
%SCALING-----------------
Ds=DiffusionCoefficient(750+273.15,XH2O);
Dscale=Ds(1);
tscale=length^2/Dscale; %dimentionless scale for the time
time=time/tscale;
nt=numel(time);
%% END:SCALING-----------------

%Initial Conditions
matrixes; % storage matrixes, separate file, KEEP IN THE SAME DIRECTORY AS THIS PROGRAM!
t=time(1)/tscale;
ZirconRadius=2e-4; %zircon radius in cm, can be changed by the user
Xs=ZirconRadius/length; %initial crystal boundary coordinate
ZircNuc=ZircNuc/length;
S=(Xs^3+MeltFrac(1)*(1-Xs^3))^(1/3);
S0=S;
dt=time(2)-time(1);
W=0; %W growth of major minerals rate
V=0; %V - zircon growth rate

C0(1:n,1)=ZrSaturation(Temp(1))*0.5145;
C0(1:n,2)=ZrSaturation(Temp(1))*0.1122;
C0(1:n,3)=ZrSaturation(Temp(1))*0.1715;
C0(1:n,4)=ZrSaturation(Temp(1))*0.1738;
C0(1:n,5)=ZrSaturation(Temp(1))*0.0280;
C0(1:n,6)=CZirc/kdHf(Temp(1),par)/70;
% C0(1:n,6)=50;  %%PHOSHPORUS< CHANGEHF melt from Bachmann etal JPet 2002.
Dplag(1:5)=0.1;
Dplag(6)=0.1;
pause(1e-5)
XXs=zeros(nt,1);
CC(1,1:n,1:5)=C0(1:n,1:5);
Tsave(1)=Temp(1)-273.15;
XXs(1)=Xs*1e4*length;
RRd(1)=S*1e4*length;
ZrPls(1)=XXs(1); %zircon radius in um
UU(1)=C0(1,1);
tt(1)=time(1)*tscale;
Zcomp(1)=C0(1,4)/C0(1,1);
ZrHF(1)=CZirc/kdHf(Temp(1),par)/C0(1,6);
Melndelta(1)=Zcomp(1);

R=linspace(0,1,n);
CZircon(1:5)=4*pi*CZirc*C0(1,1:5)/ZrSaturation(Temp(1))*ZirconRadius^3/3;
Cplag(1:5)=0;
CintS(1,1:5)=CZircon(1:5)+4*pi*C0(1,1:5)*(S^3-ZirconRadius^3)/3;
%% MAIn LOOP in time _______________________

for i=2: nt-1
    if(MeltFrac(i)>0.01)
        C=progonka(C0,dt,i,par);
        dt=time(i)-time(i-1);
        C0=C;

    else
        V=0;
        W=0;
    end
    rr=R*(S-Xs)+Xs;
    Csat=ZrSaturation(Temp(i));
    CZircon(1:5)=CZircon(1:5)-CZirc*C(1,1:5)/Csat*4*pi*Xs^2*V*dt;
    Cplag(1:5)=Cplag(1:5)-C(end,1:5).*Dplag(1:5)*4*pi*S^2*W*dt;
    Cint(1:5)=0;
    for ik=2:n
        Cint(1:5)=Cint(1:5)+(C(ik-1,1:5)*rr(ik-1)^2+C(ik,1:5)*rr(ik)^2)/2*(rr(ik)-rr(ik-1));
    end
    
    
    Cint=4*pi*Cint+par.RhoZrM*CZircon+Cplag; 
    CintS(i,1:5)=Cint(1:5);

    if(iplot==1 && mod(i,floor(nt/10))==0 )
        figure(4)
        nexttile(1)
        rr=R*(S-Xs)+Xs;
%          plot(rr*length*1e4,sum(C(:,1:5),2)-Csat,'-') %,rr*length*1e4,C(:,4),'--',LineWidth=2)
         plot(rr*length*1e4,(C(:,4)*0.5145./C(:,1)/0.1738-1)*1000,'-',LineWidth=1.5) %,rr*length*1e4,C(:,4),'--',LineWidth=2)
        hold on
        set(gca,"LineWidth",1.5, "FontSize", 14, "FontWeight", 'bold')
        ylabel('\delta^{94/90}Zr, ‰')
        xlabel('distance ,\mum')
        nexttile(2)
        
        plot(rr*length*1e4,sum(C(:,1:5),2) ./C(:,6),'-',LineWidth=1.5) %,rr*length*1e4,C(:,4),'--',LineWidth=2)
        hold on
        set(gca,"LineWidth",1.5, "FontSize", 14, "FontWeight", 'bold')
        ylabel('Zr/Hf')
        xlabel('distance ,\mum')
        
        pause(1e-10);



    end
    t=t+dt;
    rr=R*(S-Xs)+Xs;
    Cl=trapz(rr,rr.^2.*C(:,1));
    Ch=trapz(rr,rr.^2.*C(:,4));
    Xs=max(ZircNuc,Xs-V*dt);
    S0=S;
    XXs(i)=Xs*1e4*length; %zircon radius in um
    RRd(i)=S*1e4*length; %melt cell radius in um
    VV(i)=-V*length*1e4/tscale; %*length/tscale % array of dissolution rate
    tt(i)=time(i)*tscale;
    UU(i)=C(1)-ZrSaturation(Temp(i)); %-C(n
    Tsave(i)=Temp(i)-273;
    ZrPls(i)=min(XXs(1:i,1));
    Zcompl(i)=Czl/CZirc;  %    (1-(Cz/CZirc))/(Cz/CZirc);
    Zcomph(i)=Czh/CZirc;  %    (1-(Cz/CZirc))/(Cz/CZirc);
    Zcomp(i)=C(1,4)/C(1,1);
    Melndelta(i)=Ch/Cl;
    ZrHF(i)=CZirc/kdHf(Temp(i),par)/C0(1,6);
%     ZrHF(i)=kdHf(Temp(i))*C0(1,6);
    CC(i,1:n,1:6)=C0(1:n,1:6);

end
%END:MAIn LOOP in time _______________________
if (iplot==1)
    figure(1)
    nexttile(1)
    plot(tt(1:i-1)/1e3, XXs(1:i-1),'-','LineWidth',1)
    box on
    xlabel('time, ka')
    ylabel('Zircon radius ,\mum')
    set(gca,"LineWidth",1.5, "FontSize", 14, "FontWeight", 'bold')
    text(0.025,0.06, ' (a)','FontSize',14,"FontWeight", 'bold','units','normalized')
    hold on
    nexttile(2)
    plot(tt(1:i-1)/1e3, Tsave(1:i-1),'-','LineWidth',1)
    box on
    xlabel('time, ka')
    ylabel('Temperature, ^oC')
    set(gca,"LineWidth",1.5, "FontSize", 14, "FontWeight", 'bold')
    text(0.025,0.06, ' (b)','FontSize',14,"FontWeight", 'bold','units','normalized')
    hold on
    nexttile(3)
    %                 plot(tt(1:i-1)/1e3, MeltFrac(1:i-1),'-','LineWidth',1)
    plot(XXs(1:i-1), VV(1:i-1),'-','LineWidth',1) %tt(1:i-1)/1e3
    box on
    xlabel('distance, \mum')
    ylabel('Growth rate, \mum/a')
    set(gca,"LineWidth",1.5, "FontSize", 14, "FontWeight", 'bold')
    text(0.025,0.06, ' (c)','FontSize',14,"FontWeight", 'bold','units','normalized')
    hold on
    text(0.025,0.06, ' (c)','FontSize',14,"FontWeight", 'bold','units','normalized')
    nexttile(4)
    %                 plot(tt(1:i-1)/1e3, MeltFrac(1:i-1),'-','LineWidth',1)
    DelZr(1)=0;
    DelZr(2:i-1)=(Zcomp(2:i-1)./Zcomp(2)-1)*1000;
    DelMlt=(Melndelta(2:i-1)./Melndelta(2)-1)*1000;
    plot(XXs(1:i-1),DelZr(1:i-1),'-b','LineWidth',2); %(Zcomp(2:i-1)/Zcomp(2)-1)*1e3
    box on
%     ylim([-5 5])
%     hold on
%     plot(XXs(2:i-1), DelMlt,'-r','LineWidth',2);
    xlabel('distance, \mum')
    ylabel('\delta^{94/90}Zr')
%     legend('crystal','melt')
    set(gca,"LineWidth",1.5, "FontSize", 14, "FontWeight", 'bold')
    text(0.025,0.06, ' (d)','FontSize',14,"FontWeight", 'bold','units','normalized')
    hold on
    text(0.025,0.06, ' (d)','FontSize',14,"FontWeight", 'bold','units','normalized')

    nexttile(5)
    plot(XXs(2:i-1),ZrHF(2:i-1),'-','LineWidth',1); %(Zcomp(2:i-1)/Zcomp(2)-1)*1e3
    box on

    xlabel('distance, \mum')
    ylabel('Zr/Hf')
    set(gca,"LineWidth",1.5, "FontSize", 14, "FontWeight", 'bold')
    text(0.025,0.06, ' (f)','FontSize',14,"FontWeight", 'bold','units','normalized')
    hold on
    text(0.025,0.06, ' (f)','FontSize',14,"FontWeight", 'bold','units','normalized')

    nexttile(6)
    %     plot(XXs(2:i-1),(Melndelta(2:i-1)./Melndelta(2)-1)*1000 ,'-','LineWidth',1); %(Zcomp(2:i-1)/Zcomp(2)-1)*1e3
    par.Tsat=par.Tsat-273.15;
    text(0.025,0.4, formattedDisplayText(par,LineSpacing="compact"),'FontSize',10,"FontWeight", 'normal','units','normalized')
    set(gca,'visible','off')

    drawnow;
end
set(fig_th,'PaperSize',[20 10]); %set the paper size to what you want
print(fig_th,[currentFolder,tmpName,'.pdf'],'-dpdf','-fillpage') % then print it
i=nt-2;
Rsave = array2table([tt(1:i-1)/1e3,XXs(1:i-1),VV(1:i-1),Tsave(1:i-1),DelZr(1:i-1)',DelMlt(1:i-1)',ZrHF(1:i-1)']);
Rsave.Properties.VariableNames(1:7) = {'time_ka','Rad_um','Gr_rate_mm_a','Temp_C','DelZr','DelMlt','ZrHf'};
writetable(Rsave,[currentFolder,tmpName,'.csv'])
par.fname=tmpName;
if exist ('Results/summary.csv', 'file')==0
    wwar=true;
else
    wwar=false;
end

writetable(struct2table(par), 'Results/summary.csv','WriteMode','Append','Delimiter',',',...
    'WriteVariableNames',wwar,'WriteRowNames',true)
