%Script that plots (the most relevant) results

%Nominal results from Santander et al. (2022)
clear, close all
load('cases\FCC\nominalFCC_Frac_SPtrack_nom.mat');
nom.out1=R7;                                                                                                                                                                                                                     
outss=SPfcc;
nom.out2=R10;
nom.out3=R6;
nom.out4=R9;
SPT=SPfrac(:,2);
SPH1=SPfrac(:,1);
SPC1=SPfrac(:,3); %HN
SPC2=SPfrac(:,4); %LCO
out2ss=[SPH1,SPT,SPC1,SPC2];
nom.out5=[Holdup(1,:)',(Temperature(2,:)'-273.15)*(9/5)+32,(Ttrack(1,:)'),(Ttrack(2,:)')];
nom.out6=[ProductsMB(end,:)',ProductsMB(4,:)',ProductsMB(5,:)',ProductsMB(6,:)'];
nom.out7=Valvesfrac;
save('cases\FCC\results\data\FCC_Frac_SPtrack_comp.mat','nom','outss','out2ss','time');

%Nominal results from LMI tuning (This Work)
clear
load('cases\FCC\results\data\FCC_Frac_SPtrack_lmi.mat');
tw.out1=R7;                                                                                                                                                                                                                     
tw.out2=R10;
tw.out3=R6;
tw.out4=R9;
tw.out5=[Holdup(1,:)',(Temperature(2,:)'-273.15)*(9/5)+32,(Ttrack(1,:)'),(Ttrack(2,:)')];
tw.out6=[ProductsMB(end,:)',ProductsMB(4,:)',ProductsMB(5,:)',ProductsMB(6,:)'];
tw.out7=Valvesfrac;
save('cases\FCC\results\data\FCC_Frac_SPtrack_comp.mat','tw','-append');

%Nominal results from iterative LMI tuning (ILMI)
clear
load('cases\FCC\results\data\FCC_Frac_SPtrack_ilmi.mat');
ilmi.out1=R7;                                                                                                                                                                                                                     
ilmi.out2=R10;
ilmi.out3=R6;
ilmi.out4=R9;
ilmi.out5=[Holdup(1,:)',(Temperature(2,:)'-273.15)*(9/5)+32,(Ttrack(1,:)'),(Ttrack(2,:)')];
ilmi.out6=[ProductsMB(end,:)',ProductsMB(4,:)',ProductsMB(5,:)',ProductsMB(6,:)'];
ilmi.out7=Valvesfrac;
save('cases\FCC\results\data\FCC_Frac_SPtrack_comp.mat','ilmi','-append');

clear
load('cases\FCC\results\data\FCC_Frac_SPtrack_comp.mat');

%Axis labels
zfcc={'Tpre (F)';'Pfra (psig)';'Preg (psig)';'Trea (F)';'Treg (F)';'Lrea (lb)';};
zfcc_name={...
    '(a) Preheater output temperature';...
    '(d) Fractionator overhead pressure';...
    '(b) Regenerator pressure';...
    '(e) Reactor temperature';...
    '(c) Regenerator temperature';...
    '(f) Reactor catalyst inventory';};
zfcc_tick={...
    (614:2:620)';...
    (24.6:0.2:25.2)';...
    (27.6:0.2:28.2)';...
    (966:2:972)';...
    (1249:0.5:1250.5)';...
    (9.6e4:0.2e4:10.2e4)';};

MVfcc={'Ff (scf/min)';'Frgc (lb/min)';'Fsc (lb/min)';'Fa (lb/min)';'Ffg (mol/min)';'Flpg (mol/min)';};
MVfcc_name={...
    '(a) Preheater fuel flow rate';...
    '(d) Regenerated catalyst flow rate';...
    '(b) Spent catalyst flow rate';...
    '(e) Combustion air regenerator flow rate';...
    '(c) Regenerator flue gas flow rate';...
    '(f) LPG flow rate';};
MVfcc_tick={...
    (1800:200:2400)';...
    (4.9e4:0.1e4:5.2e4)';...
    (4.9e4:0.1e4:5.2e4)';...
    (3700:100:4000)';...
    (155:5:170)';...
    (2.6e4:0.1e4:2.9e4)';};

% Valvesfcc={'V4 (%)';'V6 (%)';'V7 (%)';'V3 (%)';'V1 (%)';'V2 (%)';};
% VI={'CAB (Amp)';'WGC (Amp)';};
% Prod={'VGO (lb/min)';'LPG (lb/min)';'LN (lb/min)';'HN (lb/min)';'LCO (lb/min)';'Slurry (lb/min)'};
% Valvesf={'V9 (%)';'V8 (%)';'V10 (%)';'V11 (%)'};

zf={'Lfra (kmol)';'Tfra (F)';'Thnt (F)';'Tlcot (F)';};
zf_name={...
    '(a) Fractionator accumulator level';...
    '(c) Fractionator overhead temperature';...
    '(b) HN 98% cut point';...
    '(d) LCO 98% cut point';};
zf_tick={...
    (69.8:0.1:70.1)';...
    (245.5:0.3:246.4)';...
    (527:2:533)';...
    (754:1:757)';};

MVf={'Fr (lb/min)';'Fln (lb/min)';'Fhn (lb/min)';'Flco (lb/min)';};
MVf_name={...
    '(a) Reflux flow rate';...
    '(c) LN flow rate';...
    '(b) HN flow rate';...
    '(d) LCO flow rate';};
MVf_tick={...
    (2800:100:3100)';...
    (3950:50:4100)';...
    (680:20:740)';...
    (1600:20:1660)';};

xx={'min'};
%--------------------------------------------------------------------------
%----- Plot FCC PV variables ----------------------------------------------
figure(1)
for i=1:6
    ax=subplot(3,2,i);
    if i==6
        plot(ax,...
            nom.out1(:,1),nom.out1(:,i+1)*(1/1),'-',...
            tw.out1(:,1),tw.out1(:,i+1)*(1/1),'-',...
            ilmi.out1(:,1),ilmi.out1(:,i+1)*(1/1),'-',...
            nom.out1(:,1),outss(:,i)*(1/1),'--k');
        title(zfcc_name{i},'FontName','Helvetica','FontSize',6);
        ylabel(zfcc(i)); xlabel(xx);
        set(ax,'XLim',[0 nom.out1(end,1)],'XTick',(0:50:nom.out1(end,1))','YLim',[zfcc_tick{i}(1) zfcc_tick{i}(end)],'YTick',zfcc_tick{i},'FontName','Helvetica','FontSize',6);
    else
        plot(ax,...
            nom.out1(:,1),nom.out1(:,i+1),'-',...
            tw.out1(:,1),tw.out1(:,i+1),'-',...
            ilmi.out1(:,1),ilmi.out1(:,i+1),'-',...
            nom.out1(:,1),outss(:,i),'--k');
        ylabel(zfcc(i));xlabel(xx);
        title(zfcc_name{i},'FontName','Helvetica','FontSize',6);
        ylabel(zfcc(i)); xlabel(xx);
        set(ax,'XLim',[0 nom.out1(end,1)],'XTick',(0:50:nom.out1(end,1))','YLim',[zfcc_tick{i}(1) zfcc_tick{i}(end)],'YTick',zfcc_tick{i},'FontName','Helvetica','FontSize',6);
    end
end

figh=figure(1);
figh.Units='centimeters';
figh.Renderer='painters';
figh.Color='white';
figh.OuterPosition=[0 0 14 8];
figh.Position=[0 0 14 8];

h=legend('Original','TW','ILMI');
legend('boxoff');
h.Position=[0.40,0.955,0.2,0.05];
h.Orientation='horizontal';
h.FontName='Helvetica';
h.FontSize=6;

print(figh,'figSP_PV_FCC','-dsvg');
print(figh,'figSP_PV_FCC','-depsc');

% print(figh,'figDR_PV_FCC','-dsvg');
% print(figh,'figDR_PV_FCC','-depsc');
%--------------------------------------------------------------------------
%----- Plot FCC MV variables ----------------------------------------------
figure(2)
for i=1:6
    ax=subplot(3,2,i);
    if i==6
        plot(ax,...
            nom.out2(:,1),nom.out2(:,i+1)*(1/1),'-',...
            tw.out2(:,1),tw.out2(:,i+1)*(1/1),'-',...
            ilmi.out2(:,1),ilmi.out2(:,i+1)*(1/1),'-');
        title(MVfcc_name{i},'FontName','Helvetica','FontSize',6);
        ylabel(MVfcc(i)); xlabel(xx);
        set(ax,'XLim',[0 nom.out2(end,1)],'XTick',(0:50:nom.out2(end,1))','YLim',[MVfcc_tick{i}(1) MVfcc_tick{i}(end)],'YTick',MVfcc_tick{i},'FontName','Helvetica','FontSize',6);
    else
        plot(ax,...
            nom.out2(:,1),nom.out2(:,i+1),'-',...
            tw.out2(:,1),tw.out2(:,i+1),'-',...
            ilmi.out2(:,1),ilmi.out2(:,i+1),'-');
        ylabel(MVfcc(i)); xlabel(xx);
        title(MVfcc_name{i},'FontName','Helvetica','FontSize',6);
        ylabel(MVfcc(i)); xlabel(xx);
        set(ax,'XLim',[0 nom.out2(end,1)],'XTick',(0:50:nom.out2(end,1))','YLim',[MVfcc_tick{i}(1) MVfcc_tick{i}(end)],'YTick',MVfcc_tick{i},'FontName','Helvetica','FontSize',6);
    end
end

figh=figure(2);
figh.Units='centimeters';
figh.Renderer='painters';
figh.Color='white';
figh.OuterPosition=[0 0 14 8];
figh.Position=[0 0 14 8];

h=legend('Original','TW','ILMI');
legend('boxoff');
h.Position=[0.40,0.955,0.2,0.05];
h.Orientation='horizontal';
h.FontName='Helvetica';
h.FontSize=6;

print(figh,'figSP_MV_FCC','-dsvg');
print(figh,'figSP_MV_FCC','-depsc');

% print(figh,'figDR_MV_FCC','-dsvg');
% print(figh,'figDR_MV_FCC','-depsc');
%--------------------------------------------------------------------------
%----- Plot Fractionator PV variables -------------------------------------
figure(3)
for i=1:4
    ax=subplot(2,2,i);
    plot(ax,...
        time,nom.out5(:,i),'-',...
        time,tw.out5(:,i),'-',...
        time,ilmi.out5(:,i),'-',...
        time,out2ss(:,i),'--k');
    ylabel(zf(i));xlabel(xx);
    title(zf_name{i},'FontName','Helvetica','FontSize',6);
    ylabel(zf(i)); xlabel(xx);
    set(ax,'XLim',[0 time(end,1)],'XTick',(0:50:time(end,1))','YLim',[zf_tick{i}(1) zf_tick{i}(end)],'YTick',zf_tick{i},'FontName','Helvetica','FontSize',6);
end

figh=figure(3);
figh.Units='centimeters';
figh.Renderer='painters';
figh.Color='white';
figh.OuterPosition=[0 0 14 8];
figh.Position=[0 0 14 8];

h=legend('Original','TW','ILMI');
legend('boxoff');
h.Position=[0.40,0.955,0.2,0.05];
h.Orientation='horizontal';
h.FontName='Helvetica';
h.FontSize=6;

print(figh,'figSP_PV_fra','-dsvg');
print(figh,'figSP_PV_fra','-depsc');

% print(figh,'figDR_PV_fra','-dsvg');
% print(figh,'figDR_PV_fra','-depsc');
%--------------------------------------------------------------------------
%----- Plot Fractionator MV variables -------------------------------------
figure(4)
for i=1:4
    ax=subplot(2,2,i);
    plot(ax,...
        time,nom.out6(:,i),'-',...
        time,tw.out6(:,i),'-',...
        time,ilmi.out6(:,i),'-');
    ylabel(MVf(i));xlabel(xx);
    title(MVf_name{i},'FontName','Helvetica','FontSize',6);
    ylabel(MVf(i)); xlabel(xx);
    set(ax,'XLim',[0 time(end,1)],'XTick',(0:50:time(end,1))','YLim',[MVf_tick{i}(1) MVf_tick{i}(end)],'YTick',MVf_tick{i},'FontName','Helvetica','FontSize',6);
end

figh=figure(4);
figh.Units='centimeters';
figh.Renderer='painters';
figh.Color='white';
figh.OuterPosition=[0 0 14 8];
figh.Position=[0 0 14 8];

h=legend('Original','TW','ILMI');
legend('boxoff');
h.Position=[0.40,0.955,0.2,0.05];
h.Orientation='horizontal';
h.FontName='Helvetica';
h.FontSize=6;

print(figh,'figSP_MV_fra','-dsvg');
print(figh,'figSP_MV_fra','-depsc');

% print(figh,'figDR_MV_fra','-dsvg');
% print(figh,'figDR_MV_fra','-depsc');
