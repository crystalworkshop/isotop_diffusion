function C=progonka(C0,dt,it,par)
global  n R  Dplag ZrPl MinCore time tscale S0
global A B D F alpha1 beta Xs Temp MeltFrac XH2O Tsolidus V W Csupsat Dscale UCR CZirc S ZircNuc Czl Czh

S=(Xs^3+MeltFrac(it)*(1-Xs^3))^(1/3); % rad of the melt shell
Dif=DiffusionCoefficient(Temp(it),XH2O)/Dscale; %see below Diff Coeff dependednt on water and T in cm2/s
Csat=ZrSaturation(Temp(it));
Czl=CZirc*C0(1,1)/Csat;
Czh=CZirc*C0(1,4)/Csat;

Dflux(1:5)=Dif(1:5).*(C0(2,1:5) - C0(1,1:5))/(R(2)-R(1))/(S-Xs);
V=-sum(Dflux)/(CZirc*par.RhoZrM-Csat);

if(it>1)
    diffF=(MeltFrac(it+1)-MeltFrac(it-1))/dt/2;
else
    diffF=(MeltFrac(it+1)-MeltFrac(it))/dt;
end

W=(1/3)*(diffF*(1-Xs^3)-3*Xs^2*V*(MeltFrac(it)-1))/((-MeltFrac(it)+1)*Xs^3+MeltFrac(it))^(2/3);
dC=sum(C0(n,1:5))-Csat;
Ccr=par.Crit;
delta=par.delta;
Dpmax=par.Kmax;
Dpmin=par.Kmin;
t4 = tanh(delta * (dC - Ccr));
t7 = tanh(delta * Ccr);
Dplag(1:5) = 0.1e1 / (0.1e1 + t7) * (t4 * (Dpmax - Dpmin) + Dpmax * t7 + Dpmin);
Dplag(6)=par.Ktrace;



D(n,:)=-Dif(:)-W*(R(n)-R(n-1))*(S-Xs)*(1-Dplag(:));
A(n,:)=Dif(:);
F(n,:)=0;

%Coefficients for Thomas method
s=Xs;
for j=1:6
    for i=2:n-1
        psi1=R(i-1);
        psi2=R(i);
        psi3=R(i+1);
        t1 = (Dif(j) * dt);
        t5 = (psi1 * S - psi1 * s + s) ^ 2;
        t6 = psi2 - psi1;
        t8 = (t5 / t6);
        t12 = (S * psi2);
        t14 = ((-psi2 + 1) * s + t12) ^ 2;
        t15 = S - s;
        t20 = (-W + V) * psi2 - V;
        A(i,j) = -t14 * t15 * dt * psi2 * t20 - t1 * t8;
        t25 = (-psi2 * s + s + t12) ^ 2;
        t28 = t25 / (psi3 - psi2);
        B(i,j) = -t1 * t28;
        t32 = -t15;
        t33 = t32 ^ 2;
        t34 = -t6;
        t38 = (t32 * psi2 - s) ^ 2;
        D(i,j) = -t1 * (-t28 - t8) - t33 * t34 * t38 - t20 * psi2 * dt * t38 * t32;
        t44 = t15 ^ 2;
        t48 = (t15 * psi2 + s) ^ 2;
        F(i,j) = -t34 * t44 * t48 * C0(i,j);
    end
end
%Forward Thomas path
alpha1(n,:)=-A(n,:)./D(n,:);
beta(n,:)=F(n,:)/D(n,:);
for i=n-1:-1:2
    alpha1(i,:)=-A(i,:)./(B(i,:).*alpha1(i+1,:)+D(i,:));
    beta(i,:)=(F(i,:)-B(i,:).*beta(i+1,:))./(B(i,:).*alpha1(i+1,:)+D(i,:));
end
% Boundary conditions
parb.D=Dif(:);
parb.csat=Csat;
parb.alpha=alpha1(2,:);
parb.beta=beta(2,:);
parb.T=Temp(it);
parb.Cz=CZirc;
parb.Trace=par.Trace;
f = @(X) bc(X,parb); % function of dummy variable y
[out,fval,exflag]=fsolve(f,C0(1,:),optimset('Display','off','Algorithm','levenberg-marquardt'));
if exflag<=0 
    disp(fval)
end
C(1,:)=out(:);


%Backward Thomas path
for i=1:n-1
    C(i+1,:)=C(i,:).*alpha1(i+1,:)+beta(i+1,:);
end

end



