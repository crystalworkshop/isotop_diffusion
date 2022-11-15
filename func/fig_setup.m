
f = figure('Units','normalized','Visible','off','Position',[.10,.10,.8,.8],'Color',[1 1 1]);
 
ha = axes('Units','normalized','Position',[0.05,.38,.40,0.57]);
xlim([0 1]); box on
ha1 = axes('Units','normalized','Position',[0.05,0.05,0.40,0.28]);
%xlim([0 1]); 
box on
 
ha2 = axes('Units','normalized','Position',[0.5,0.65,0.40,0.30]); box on
ha3 = axes('Units','normalized','Position',[0.5,0.45,0.40,0.15]); box on
ha4 = axes('Units','normalized','Position',[0.5,0.25,0.40,0.15]); box on
ha5 = axes('Units','normalized','Position',[0.5,0.05,0.40,0.15]); box on
 
xlabel(ha1,'Distance, um')
ylabel(ha1,'Zr concentration')
xlabel(ha2,'Time (years)')
ylabel(ha2,'Zr radius')
xlabel(ha3,'Time (years)')
ylabel(ha3,'Growth Rate, cm.s^{-1}')
xlabel(ha4,'Time (years)')
ylabel(ha4,'Temperature T, ^oC')
xlabel(ha5,'Time (years)')
ylabel(ha5,'Crystal content')
 
 
set(f,'CurrentAxes',ha)
set(f,'Name','DIFFUSOR')
movegui(f,'center')
set(f,'Visible','on');
 
hold on
tyear=365*3600*24;
