function [ti, Ti,CrFrac1]=TemperatureHistory_m_erupt(tr,Tr,par)
if isempty(tr)
    nt=par.nt;
    ti=linspace(0,par.tfin,nt);
    Ti=linspace(par.Tsat, par.Tend+273.15,nt);
    CrFrac1=mf_rock(Ti-273.15); 
   % CrFrac1=linspace(1,1,nt); %mf_rock(Ti-273.15); 
else
istart=find(Tr>0,1);
Tr(istart-1)=950+273.15;
tr(istart-1)=tr(istart)-5;
dT=0.05;
if min(mf_rock(Tr))<0.01
    disp('no melting')
    ti=[]; Ti=[]; CrFrac1=[];
    return
end
if  min(Tr)-Tsat >0
    disp('high temperature')
    ti=[]; Ti=[]; CrFrac1=[];
    return
end
time=tr;
Temp=Tr;
try
it=find(Temp<Tsat,1);
time(it-1)=time(it)-(Temp(it)-Tsat)/(Temp(it)-Temp(it-1))*5;
Temp(it-1)=Tsat;
time1=time(it-1:end);
Temp1=Temp(it-1:end);
nt=numel(time1);
s=zeros(size(Temp1));
for i=2:nt
    s(i)=s(i-1)+abs(Temp1(i)-Temp1(i-1));
end

ni=floor(s(nt)/dT);
si=linspace(s(1),s(nt),ni);
ti=interp1(s,time1,si);
Ti=interp1(time1,Temp1,ti);
catch ME
    warning(['wrong Thist for sample: ',num2str(sampnum), ',',ME.message]);
    ti=[]; Ti=[]; CrFrac1=[]; return
end
CrFrac1=mf_rock(Ti-273.15);
end
end
