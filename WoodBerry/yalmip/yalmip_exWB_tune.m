clear, close all
load('cases\WoodBerry\yalmip\results\data\yalmip_exWB_tb_rlt.mat');
load(rlt(1).data,'A_pi', 'Bu_pi', 'Bw_pi', 'Cy_pi', 'A_pid', 'Bu_pid', 'Bw_pid', 'Cy_pid', 'Bu', 'Bw', 'Cy', 'PhiD', 'Nx', 'Ny', 'Nw');

tstart=0; tend=300; tstep=1; t=(tstart:tstep:tend)'; %time domain
yf=[0.1;0.1]; dw=0.245;

%Tuning #1 from this work
KP_TW1=rlt(1).TW.PI.KP;
KI_TW1=rlt(1).TW.PI.KI;

%Tuning #5 from this work
KP_TW5=rlt(3).TW.PI.KP;
KI_TW5=rlt(3).TW.PI.KI;

%Tuning #7 from this work
KP_TW7=rlt(1).TW.PID.KP;
KI_TW7=rlt(1).TW.PID.KI;
KD_TW7=rlt(1).TW.PID.KD;

%Tuning #12 from this work
KP_TW12=rlt(3).TW.PID.KP;
KI_TW12=rlt(3).TW.PID.KI;
KD_TW12=rlt(3).TW.PID.KD;

%Tuning from Boyd, S., Hast, M., Astrom, J. (2016)
KP_BHA16=...
	[0.1535 	  0
	 0      -0.0692];
	 
KI_BHA16=...
	[0.0210 	  0
	 0      -0.0136];

KD_BHA16=...
	[0.1714 	  0
	 0      -0.1725];
	 
%Tuning from Chen, D., Seaborg, D. E. (2003)
KP_CS03=...
	[0.436 	  0
	 0      -0.0945];
	 
KI_CS03=...
	[11\0.436 	  0
	 0      -15.5\0.0945];
	 
%BLT tuning Luyben, W. L. (1986)
KP_BLT=...
	[0.375 	  0
	 0      -0.075];
	 
KI_BLT=...
	[8.29\0.375 	  0
	 0      -23.6\0.075];	 

%ITAE based tuning for FOPDT
KP_ITAE=[
    (1/12.8)*0.965*(1/16.7)^(-0.855);...
    (1/-19.4)*0.965*(3/14.4)^(-0.855)];
KP_ITAE=diag(KP_ITAE);

TI_ITAE=[
    16.7/(0.796+0.147*(1/16.7));...
    14.4/(0.796+0.147*(3/14.4))];
KI_ITAE=diag(TI_ITAE)\KP_ITAE;

TD_ITAE=[
    16.7*(0.308*(1/16.7)^0.9292);...
    14.4*(0.308*(3/14.4)^0.9292)];
KD_ITAE=KP_ITAE*diag(TD_ITAE);


%Tuning from Garrido et al. (2021): https://doi.org/10.3390/pr9010140
KP_GAR21=...
	[0.6500 	  0
	 0      -0.0930];
	 
KI_GAR21=...
	[0.0640 	  0
	 0      -0.0320];

KD_GAR21=...
	[0.2600 	  0
	 0      -0.1760];

%----- Plot PV dynamic response -------------------------------------
%Plot step response
dataopt=stepDataOptions('StepAmplitude',[yf;dw]);

%This work #1
Br_=[Bu*KP_TW1;
	 -eye(Ny)];
Bw_=[Bw
	 zeros(Ny,Nw)];
sys_TW1=ss(A_pi-Bu_pi*[KP_TW1,KI_TW1]*Cy_pi, [Br_,Bw_], Cy_pi, []);	 
step_TW1=step(sys_TW1,t,dataopt);
yPV_TW1(:,1)=step_TW1(:,1,1);
yPV_TW1(:,2)=step_TW1(:,1,2);
yPV_TW1(:,3)=step_TW1(:,1,3);
yPV_TW1(:,4)=step_TW1(:,2,1);
yPV_TW1(:,5)=step_TW1(:,2,2);
yPV_TW1(:,6)=step_TW1(:,2,3);

%This work #5
Br_=[Bu*KP_TW5;
	 -eye(Ny)];
Bw_=[Bw
	 zeros(Ny,Nw)];
sys_TW5=ss(A_pi-Bu_pi*[KP_TW5,KI_TW5]*Cy_pi, [Br_,Bw_], Cy_pi, []);	 
step_TW5=step(sys_TW5,t,dataopt);
yPV_TW5(:,1)=step_TW5(:,1,1);
yPV_TW5(:,2)=step_TW5(:,1,2);
yPV_TW5(:,3)=step_TW5(:,1,3);
yPV_TW5(:,4)=step_TW5(:,2,1);
yPV_TW5(:,5)=step_TW5(:,2,2);
yPV_TW5(:,6)=step_TW5(:,2,3);

%This work #7
Br_=[Bu*KP_TW7;
	 -eye(Ny);
	 PhiD*Cy*Bu*KP_TW7];
Bw_=[Bw
	 zeros(Ny,Nw)
	 PhiD*Cy*Bw];
sys_TW7=ss(A_pid-Bu_pid*[KP_TW7,KI_TW7,KD_TW7]*Cy_pid, [Br_,Bw_], Cy_pid, []);
step_TW7=step(sys_TW7,t,dataopt);
yPV_TW7(:,1)=step_TW7(:,1,1);
yPV_TW7(:,2)=step_TW7(:,1,2);
yPV_TW7(:,3)=step_TW7(:,1,3);
yPV_TW7(:,4)=step_TW7(:,2,1);
yPV_TW7(:,5)=step_TW7(:,2,2);
yPV_TW7(:,6)=step_TW7(:,2,3);

%This work #12
Br_=[Bu*KP_TW12;
	 -eye(Ny);
	 PhiD*Cy*Bu*KP_TW12];
Bw_=[Bw
	 zeros(Ny,Nw)
	 PhiD*Cy*Bw];
sys_TW12=ss(A_pid-Bu_pid*[KP_TW12,KI_TW12,KD_TW12]*Cy_pid, [Br_,Bw_], Cy_pid, []);
step_TW12=step(sys_TW12,t,dataopt);
yPV_TW12(:,1)=step_TW12(:,1,1);
yPV_TW12(:,2)=step_TW12(:,1,2);
yPV_TW12(:,3)=step_TW12(:,1,3);
yPV_TW12(:,4)=step_TW12(:,2,1);
yPV_TW12(:,5)=step_TW12(:,2,2);
yPV_TW12(:,6)=step_TW12(:,2,3);

%Boyd, Hast & Astrom (2016)
Br_=[Bu*KP_BHA16;
	 -eye(Ny);
	 PhiD*Cy*Bu*KP_BHA16];
Bw_=[Bw
	 zeros(Ny,Nw)
	 PhiD*Cy*Bw];
sys_BHA16=ss(A_pid-Bu_pid*[KP_BHA16,KI_BHA16,KD_BHA16]*Cy_pid, [Br_,Bw_], Cy_pid, []);
step_BHA16=step(sys_BHA16,t,dataopt);
yPV_BHA16(:,1)=step_BHA16(:,1,1);
yPV_BHA16(:,2)=step_BHA16(:,1,2);
yPV_BHA16(:,3)=step_BHA16(:,1,3);
yPV_BHA16(:,4)=step_BHA16(:,2,1);
yPV_BHA16(:,5)=step_BHA16(:,2,2);
yPV_BHA16(:,6)=step_BHA16(:,2,3);

%Chen & Seaborg (2003)
Br_=[Bu*KP_CS03;
	 -eye(Ny)];
Bw_=[Bw
	 zeros(Ny,Nw)];
sys_CS03=ss(A_pi-Bu_pi*[KP_CS03,KI_CS03]*Cy_pi,[Br_,Bw_], Cy_pi, []);
step_CS03=step(sys_CS03,t,dataopt);
yPV_CS03(:,1)=step_CS03(:,1,1);
yPV_CS03(:,2)=step_CS03(:,1,2);
yPV_CS03(:,3)=step_CS03(:,1,3);
yPV_CS03(:,4)=step_CS03(:,2,1);
yPV_CS03(:,5)=step_CS03(:,2,2);
yPV_CS03(:,6)=step_CS03(:,2,3);

%Biggest-Log-Modulos (Luyben, 1986)
Br_=[Bu*KP_BLT;
	 -eye(Ny)];
Bw_=[Bw
	 zeros(Ny,Nw)];
sys_BLT=ss(A_pi-Bu_pi*[KP_BLT,KI_BLT]*Cy_pi,[Br_,Bw_], Cy_pi, []);
step_BLT=step(sys_BLT,t,dataopt);
yPV_BLT(:,1)=step_BLT(:,1,1);
yPV_BLT(:,2)=step_BLT(:,1,2);
yPV_BLT(:,3)=step_BLT(:,1,3);
yPV_BLT(:,4)=step_BLT(:,2,1);
yPV_BLT(:,5)=step_BLT(:,2,2);
yPV_BLT(:,6)=step_BLT(:,2,3);

%ITAE based tuning for FOPDT
Br_=[Bu*KP_ITAE;
	 -eye(Ny);
	 PhiD*Cy*Bu*KP_ITAE];
Bw_=[Bw
	 zeros(Ny,Nw)
	 PhiD*Cy*Bw];
sys_ITAE=ss(A_pid-Bu_pid*[KP_ITAE,KI_ITAE,KD_ITAE]*Cy_pid, [Br_,Bw_], Cy_pid, []);
step_ITAE=step(sys_ITAE,t,dataopt);
yPV_ITAE(:,1)=step_ITAE(:,1,1);
yPV_ITAE(:,2)=step_ITAE(:,1,2);
yPV_ITAE(:,3)=step_ITAE(:,1,3);
yPV_ITAE(:,4)=step_ITAE(:,2,1);
yPV_ITAE(:,5)=step_ITAE(:,2,2);
yPV_ITAE(:,6)=step_ITAE(:,2,3);

%Garrido et al. (2021)
Br_=[Bu*KP_GAR21;
	 -eye(Ny);
	 PhiD*Cy*Bu*KP_GAR21];
Bw_=[Bw
	 zeros(Ny,Nw)
	 PhiD*Cy*Bw];
sys_GAR21=ss(A_pid-Bu_pid*[KP_GAR21,KI_GAR21,KD_GAR21]*Cy_pid, [Br_,Bw_], Cy_pid, []);
step_GAR21=step(sys_GAR21,t,dataopt);
yPV_GAR21(:,1)=step_GAR21(:,1,1);
yPV_GAR21(:,2)=step_GAR21(:,1,2);
yPV_GAR21(:,3)=step_GAR21(:,1,3);
yPV_GAR21(:,4)=step_GAR21(:,2,1);
yPV_GAR21(:,5)=step_GAR21(:,2,2);
yPV_GAR21(:,6)=step_GAR21(:,2,3);

%Axis labels
out={'Distillate [wt. frac.]';'';'';'Bottom [wt. frac.]';'';''};
out_name={'Dist. wt. frac. SP';...
    'Bottom wt. frac. SP';...
	'Feed';...
    '';...
    '';...
    ''};
out_tick={...
    (-0.2:0.08:0.12)';...
    (-0.2:0.08:0.12)';...
    (-0.2:0.08:0.12)';...
    (-0.2:0.8/4:0.6)';...
    (-0.2:0.8/4:0.6)';...
    (-0.2:0.8/4:0.6)'};
xx={'min'};

figh=figure(1);
for i=1:6
ax=subplot(2,3,i);
plot(ax,...
    ...t,yPV_TW1(:,i),'-',...
    t,yPV_TW5(:,i),'-',...
    ...t,yPV_TW7(:,i),'-',...
    ...t,yPV_TW12(:,i),'-',...
    t,yPV_BHA16(:,i),'-',...
	...t,yPV_CS03(:,i),'-',...
	...t,yPV_BLT(:,i),'-',...
    t,yPV_ITAE(:,i),'-',...
    t,yPV_GAR21(:,i),'-',...
    t,ones(size(t,1))*yPV_TW1(end,i),'--k');
title(out_name{i},'FontName','Helvetica','FontSize',6);
ylabel(out(i),'FontName','Helvetica','FontSize',6);
xlabel(xx);
set(ax,'XLim',[0 t(end,1)],'XTick',(0:60:t(end,1))','YLim',[out_tick{i}(1) out_tick{i}(end)],'YTick',out_tick{i},'FontName','Helvetica','FontSize',6);
end

figh.Units='centimeters';
figh.Color='none'; % figh.Color='white';
figh.Renderer='painters';
figh.OuterPosition=[0 0 14 8];
figh.Position=[0 0 14 8];

% h=legend('TW1','TW5','TW7','BHA16','CS03','BLT');
% h=legend('TW','BHA16','CS03','BLT');
% h=legend('TW','BHA16','CS03','BLT','ITAE-SP');
% h=legend('TW','BHA16','CS03','ITAE-SP');
% h=legend('TW','BHA16','CS03','ITAE-SP','GAR21');
h=legend('TW','BHA16','ITAE-SP','GAR21');
legend('boxoff');
h.Position=[0.40,0.955,0.2,0.05];
h.Orientation='horizontal';
h.FontName='Helvetica';
h.FontSize=6;

print(figh,'fig_WBtune_PV','-dsvg');
print(figh,'fig_WBtune_PV','-depsc');
%--------------------------------------------------------------------
%----- Plot MV dynamic response -------------------------------------
%Reference input
yr(:,:,1)=ones(size(t,1),1)*[yf(1);0]';
yr(:,:,2)=ones(size(t,1),1)*[0;yf(2)]';
yr(:,:,3)=ones(size(t,1),1)*[0;0]';

%This work #1
uMV_TW1(:,1)=([KP_TW1(1,1) -[KP_TW1(1,1),KI_TW1(1,1)]]*[yr(:,1,1), step_TW1(:,[1;3],1)]')';
uMV_TW1(:,2)=([KP_TW1(1,1) -[KP_TW1(1,1),KI_TW1(1,1)]]*[yr(:,1,2), step_TW1(:,[1;3],2)]')';
uMV_TW1(:,3)=([KP_TW1(1,1) -[KP_TW1(1,1),KI_TW1(1,1)]]*[yr(:,1,3), step_TW1(:,[1;3],3)]')';
uMV_TW1(:,4)=([KP_TW1(2,2) -[KP_TW1(2,2),KI_TW1(2,2)]]*[yr(:,2,1), step_TW1(:,[2;4],1)]')';
uMV_TW1(:,5)=([KP_TW1(2,2) -[KP_TW1(2,2),KI_TW1(2,2)]]*[yr(:,2,2), step_TW1(:,[2;4],2)]')';
uMV_TW1(:,6)=([KP_TW1(2,2) -[KP_TW1(2,2),KI_TW1(2,2)]]*[yr(:,2,3), step_TW1(:,[2;4],3)]')';

%This work #5
uMV_TW5(:,1)=([KP_TW5(1,1) -[KP_TW5(1,1),KI_TW5(1,1)]]*[yr(:,1,1), step_TW5(:,[1;3],1)]')';
uMV_TW5(:,2)=([KP_TW5(1,1) -[KP_TW5(1,1),KI_TW5(1,1)]]*[yr(:,1,2), step_TW5(:,[1;3],2)]')';
uMV_TW5(:,3)=([KP_TW5(1,1) -[KP_TW5(1,1),KI_TW5(1,1)]]*[yr(:,1,3), step_TW5(:,[1;3],3)]')';
uMV_TW5(:,4)=([KP_TW5(2,2) -[KP_TW5(2,2),KI_TW5(2,2)]]*[yr(:,2,1), step_TW5(:,[2;4],1)]')';
uMV_TW5(:,5)=([KP_TW5(2,2) -[KP_TW5(2,2),KI_TW5(2,2)]]*[yr(:,2,2), step_TW5(:,[2;4],2)]')';
uMV_TW5(:,6)=([KP_TW5(2,2) -[KP_TW5(2,2),KI_TW5(2,2)]]*[yr(:,2,3), step_TW5(:,[2;4],3)]')';

%This work #7
uMV_TW7(:,1)=([KP_TW7(1,1) -[KP_TW7(1,1),KI_TW7(1,1),KD_TW7(1,1)]]*[yr(:,1,1), step_TW7(:,[1;3;5],1)]')';
uMV_TW7(:,2)=([KP_TW7(1,1) -[KP_TW7(1,1),KI_TW7(1,1),KD_TW7(1,1)]]*[yr(:,1,2), step_TW7(:,[1;3;5],2)]')';
uMV_TW7(:,3)=([KP_TW7(1,1) -[KP_TW7(1,1),KI_TW7(1,1),KD_TW7(1,1)]]*[yr(:,1,3), step_TW7(:,[1;3;5],3)]')';
uMV_TW7(:,4)=([KP_TW7(2,2) -[KP_TW7(2,2),KI_TW7(2,2),KD_TW7(2,2)]]*[yr(:,2,1), step_TW7(:,[2;4;6],1)]')';
uMV_TW7(:,5)=([KP_TW7(2,2) -[KP_TW7(2,2),KI_TW7(2,2),KD_TW7(2,2)]]*[yr(:,2,2), step_TW7(:,[2;4;6],2)]')';
uMV_TW7(:,6)=([KP_TW7(2,2) -[KP_TW7(2,2),KI_TW7(2,2),KD_TW7(2,2)]]*[yr(:,2,3), step_TW7(:,[2;4;6],3)]')';

%This work #12
uMV_TW12(:,1)=([KP_TW12(1,1) -[KP_TW12(1,1),KI_TW12(1,1),KD_TW12(1,1)]]*[yr(:,1,1), step_TW12(:,[1;3;5],1)]')';
uMV_TW12(:,2)=([KP_TW12(1,1) -[KP_TW12(1,1),KI_TW12(1,1),KD_TW12(1,1)]]*[yr(:,1,2), step_TW12(:,[1;3;5],2)]')';
uMV_TW12(:,3)=([KP_TW12(1,1) -[KP_TW12(1,1),KI_TW12(1,1),KD_TW12(1,1)]]*[yr(:,1,3), step_TW12(:,[1;3;5],3)]')';
uMV_TW12(:,4)=([KP_TW12(2,2) -[KP_TW12(2,2),KI_TW12(2,2),KD_TW12(2,2)]]*[yr(:,2,1), step_TW12(:,[2;4;6],1)]')';
uMV_TW12(:,5)=([KP_TW12(2,2) -[KP_TW12(2,2),KI_TW12(2,2),KD_TW12(2,2)]]*[yr(:,2,2), step_TW12(:,[2;4;6],2)]')';
uMV_TW12(:,6)=([KP_TW12(2,2) -[KP_TW12(2,2),KI_TW12(2,2),KD_TW12(2,2)]]*[yr(:,2,3), step_TW12(:,[2;4;6],3)]')';

%Boyd, Hast & Astrom (2016)
uMV_BHA16(:,1)=([KP_BHA16(1,1) -[KP_BHA16(1,1),KI_BHA16(1,1),KD_BHA16(1,1)]]*[yr(:,1,1), step_BHA16(:,[1;3;5],1)]')';
uMV_BHA16(:,2)=([KP_BHA16(1,1) -[KP_BHA16(1,1),KI_BHA16(1,1),KD_BHA16(1,1)]]*[yr(:,1,2), step_BHA16(:,[1;3;5],2)]')';
uMV_BHA16(:,3)=([KP_BHA16(1,1) -[KP_BHA16(1,1),KI_BHA16(1,1),KD_BHA16(1,1)]]*[yr(:,1,3), step_BHA16(:,[1;3;5],3)]')';
uMV_BHA16(:,4)=([KP_BHA16(2,2) -[KP_BHA16(2,2),KI_BHA16(2,2),KD_BHA16(2,2)]]*[yr(:,2,1), step_BHA16(:,[2;4;6],1)]')';
uMV_BHA16(:,5)=([KP_BHA16(2,2) -[KP_BHA16(2,2),KI_BHA16(2,2),KD_BHA16(2,2)]]*[yr(:,2,2), step_BHA16(:,[2;4;6],2)]')';
uMV_BHA16(:,6)=([KP_BHA16(2,2) -[KP_BHA16(2,2),KI_BHA16(2,2),KD_BHA16(2,2)]]*[yr(:,2,3), step_BHA16(:,[2;4;6],3)]')';

%Chen & Seaborg (2003)
uMV_CS03(:,1)=([KP_CS03(1,1) -[KP_CS03(1,1),KI_CS03(1,1)]]*[yr(:,1,1), step_CS03(:,[1;3],1)]')';
uMV_CS03(:,2)=([KP_CS03(1,1) -[KP_CS03(1,1),KI_CS03(1,1)]]*[yr(:,1,2), step_CS03(:,[1;3],2)]')';
uMV_CS03(:,3)=([KP_CS03(1,1) -[KP_CS03(1,1),KI_CS03(1,1)]]*[yr(:,1,3), step_CS03(:,[1;3],3)]')';
uMV_CS03(:,4)=([KP_CS03(2,2) -[KP_CS03(2,2),KI_CS03(2,2)]]*[yr(:,2,1), step_CS03(:,[2;4],1)]')';
uMV_CS03(:,5)=([KP_CS03(2,2) -[KP_CS03(2,2),KI_CS03(2,2)]]*[yr(:,2,2), step_CS03(:,[2;4],2)]')';
uMV_CS03(:,6)=([KP_CS03(2,2) -[KP_CS03(2,2),KI_CS03(2,2)]]*[yr(:,2,3), step_CS03(:,[2;4],3)]')';

%Biggest-Log-Modulos (Luyben, 1986)
uMV_BLT(:,1)=([KP_BLT(1,1) -[KP_BLT(1,1),KI_BLT(1,1)]]*[yr(:,1,1), step_BLT(:,[1;3],1)]')';
uMV_BLT(:,2)=([KP_BLT(1,1) -[KP_BLT(1,1),KI_BLT(1,1)]]*[yr(:,1,2), step_BLT(:,[1;3],2)]')';
uMV_BLT(:,3)=([KP_BLT(1,1) -[KP_BLT(1,1),KI_BLT(1,1)]]*[yr(:,1,3), step_BLT(:,[1;3],3)]')';
uMV_BLT(:,4)=([KP_BLT(2,2) -[KP_BLT(2,2),KI_BLT(2,2)]]*[yr(:,2,1), step_BLT(:,[2;4],1)]')';
uMV_BLT(:,5)=([KP_BLT(2,2) -[KP_BLT(2,2),KI_BLT(2,2)]]*[yr(:,2,2), step_BLT(:,[2;4],2)]')';
uMV_BLT(:,6)=([KP_BLT(2,2) -[KP_BLT(2,2),KI_BLT(2,2)]]*[yr(:,2,3), step_BLT(:,[2;4],3)]')';

%ITAE based tuning for FOPDT
uMV_ITAE(:,1)=([KP_ITAE(1,1) -[KP_ITAE(1,1),KI_ITAE(1,1),KD_ITAE(1,1)]]*[yr(:,1,1), step_ITAE(:,[1;3;5],1)]')';
uMV_ITAE(:,2)=([KP_ITAE(1,1) -[KP_ITAE(1,1),KI_ITAE(1,1),KD_ITAE(1,1)]]*[yr(:,1,2), step_ITAE(:,[1;3;5],2)]')';
uMV_ITAE(:,3)=([KP_ITAE(1,1) -[KP_ITAE(1,1),KI_ITAE(1,1),KD_ITAE(1,1)]]*[yr(:,1,3), step_ITAE(:,[1;3;5],3)]')';
uMV_ITAE(:,4)=([KP_ITAE(2,2) -[KP_ITAE(2,2),KI_ITAE(2,2),KD_ITAE(2,2)]]*[yr(:,2,1), step_ITAE(:,[2;4;6],1)]')';
uMV_ITAE(:,5)=([KP_ITAE(2,2) -[KP_ITAE(2,2),KI_ITAE(2,2),KD_ITAE(2,2)]]*[yr(:,2,2), step_ITAE(:,[2;4;6],2)]')';
uMV_ITAE(:,6)=([KP_ITAE(2,2) -[KP_ITAE(2,2),KI_ITAE(2,2),KD_ITAE(2,2)]]*[yr(:,2,3), step_ITAE(:,[2;4;6],3)]')';

%Garrido et al. (2021)
uMV_GAR21(:,1)=([KP_GAR21(1,1) -[KP_GAR21(1,1),KI_GAR21(1,1),KD_GAR21(1,1)]]*[yr(:,1,1), step_GAR21(:,[1;3;5],1)]')';
uMV_GAR21(:,2)=([KP_GAR21(1,1) -[KP_GAR21(1,1),KI_GAR21(1,1),KD_GAR21(1,1)]]*[yr(:,1,2), step_GAR21(:,[1;3;5],2)]')';
uMV_GAR21(:,3)=([KP_GAR21(1,1) -[KP_GAR21(1,1),KI_GAR21(1,1),KD_GAR21(1,1)]]*[yr(:,1,3), step_GAR21(:,[1;3;5],3)]')';
uMV_GAR21(:,4)=([KP_GAR21(2,2) -[KP_GAR21(2,2),KI_GAR21(2,2),KD_GAR21(2,2)]]*[yr(:,2,1), step_GAR21(:,[2;4;6],1)]')';
uMV_GAR21(:,5)=([KP_GAR21(2,2) -[KP_GAR21(2,2),KI_GAR21(2,2),KD_GAR21(2,2)]]*[yr(:,2,2), step_GAR21(:,[2;4;6],2)]')';
uMV_GAR21(:,6)=([KP_GAR21(2,2) -[KP_GAR21(2,2),KI_GAR21(2,2),KD_GAR21(2,2)]]*[yr(:,2,3), step_GAR21(:,[2;4;6],3)]')';

%Axis labels
out={'Reflux [lb/min]';'';'';'Vapor [lb/min]';'';''};
out_name={'Dist. wt. frac. SP';...
    'Bottom wt. frac. SP';...
	'Feed';...
    '';...
    '';...
    ''};
out_tick={...
    (-0.05:0.1/4:0.05)';...
    (-0.05:0.1/4:0.05)';...
    (-0.05:0.1/4:0.05)';...
    (-0.02:0.12/4:0.10)';...
    (-0.02:0.12/4:0.10)';...
    (-0.02:0.12/4:0.10)'};
xx={'min'};

figh=figure(2);
for i=1:6
ax=subplot(2,3,i);
plot(ax,...
    ...t,uMV_TW1(:,i),'-',...
    t,uMV_TW5(:,i),'-',...
    ...t,uMV_TW7(:,i),'-',...
    ...t,uMV_TW12(:,i),'-',...
    t,uMV_BHA16(:,i),'-',...
	...t,uMV_CS03(:,i),'-',...
	...t,uMV_BLT(:,i),'-',...
    t,uMV_ITAE(:,i),'-',...
    t,uMV_GAR21(:,i),'-',...
    t,ones(size(t,1))*uMV_TW1(end,i),'--k');
title(out_name{i},'FontName','Helvetica','FontSize',6);
ylabel(out(i),'FontName','Helvetica','FontSize',6);
xlabel(xx);
set(ax,'XLim',[0 t(end,1)],'XTick',(0:60:t(end,1))','YLim',[out_tick{i}(1) out_tick{i}(end)],'YTick',out_tick{i},'FontName','Helvetica','FontSize',6);
end

figh.Units='centimeters';
figh.Color='none'; % figh.Color='white'; 
figh.Renderer='painters';
figh.OuterPosition=[0 0 14 8];
figh.Position=[0 0 14 8];

% h=legend('TW1','TW5','TW7','BHA16','CS03','BLT');
% h=legend('TW','BHA16','CS03','BLT');
% h=legend('TW','BHA16','CS03','BLT','ITAE-SP');
% h=legend('TW','BHA16','CS03','ITAE-SP');
% h=legend('TW','BHA16','CS03','ITAE-SP','GAR21');
h=legend('TW','BHA16','ITAE-SP','GAR21');
legend('boxoff');
h.Position=[0.40,0.955,0.2,0.05];
h.Orientation='horizontal';
h.FontName='Helvetica';
h.FontSize=6;

print(figh,'fig_WBtune_MV','-dsvg');
print(figh,'fig_WBtune_MV','-depsc');
%--------------------------------------------------------------------
%----- Table 4 ------------------------------------------------------
%ITAE comparison
info_TW1=stepinfo_caseWB(KP_TW1,KI_TW1,[],[],ss(A_pi,[Bu_pi,Bw_pi],Cy_pi,0),t,[yf;dw],'PI');
info_TW5=stepinfo_caseWB(KP_TW5,KI_TW5,[],[],ss(A_pi,[Bu_pi,Bw_pi],Cy_pi,0),t,[yf;dw],'PI');
info_TW7=stepinfo_caseWB(KP_TW7,KI_TW7,KD_TW7,PhiD,ss(A_pid,[Bu_pid,Bw_pid],Cy_pid,0),t,[yf;dw],'PID');
info_TW12=stepinfo_caseWB(KP_TW12,KI_TW12,KD_TW12,PhiD,ss(A_pid,[Bu_pid,Bw_pid],Cy_pid,0),t,[yf;dw],'PID');
info_BHA16=stepinfo_caseWB(KP_BHA16,KI_BHA16,KD_BHA16,PhiD,ss(A_pid,[Bu_pid,Bw_pid],Cy_pid,0),t,[yf;dw],'PID');
info_CS03=stepinfo_caseWB(KP_CS03,KI_CS03,[],[],ss(A_pi,[Bu_pi,Bw_pi],Cy_pi,0),t,[yf;dw],'PI');
info_BLT=stepinfo_caseWB(KP_BLT,KI_BLT,[],[],ss(A_pi,[Bu_pi,Bw_pi],Cy_pi,0),t,[yf;dw],'PI');
info_ITAE=stepinfo_caseWB(KP_ITAE,KI_ITAE,KD_ITAE,PhiD,ss(A_pid,[Bu_pid,Bw_pid],Cy_pid,0),t,[yf;dw],'PID');
info_GAR21=stepinfo_caseWB(KP_GAR21,KI_GAR21,KD_GAR21,PhiD,ss(A_pid,[Bu_pid,Bw_pid],Cy_pid,0),t,[yf;dw],'PID');

%Disp table 5
tab5=table(...
    num2str([KP_TW1(1,1);  		            KP_TW1(2,2);...
             KP_TW5(1,1);  		            KP_TW5(2,2);...
             KP_TW7(1,1);  		            KP_TW7(2,2);...
             KP_TW12(1,1);  		        KP_TW12(2,2);...
			 KP_BHA16(1,1);  				KP_BHA16(2,2);...
			 KP_CS03(1,1);  				KP_CS03(2,2);...
			 KP_BLT(1,1);  					KP_BLT(2,2);...
             KP_ITAE(1,1);  				KP_ITAE(2,2);...
             KP_GAR21(1,1);  				KP_GAR21(2,2)],'%.4f'),...
    num2str([KI_TW1(1,1)\KP_TW1(1,1);  		KI_TW1(2,2)\KP_TW1(2,2);...
             KI_TW5(1,1)\KP_TW5(1,1);  		KI_TW5(2,2)\KP_TW5(2,2);...
             KI_TW7(1,1)\KP_TW7(1,1);  		KI_TW7(2,2)\KP_TW7(2,2);...
             KI_TW12(1,1)\KP_TW12(1,1);  	KI_TW12(2,2)\KP_TW12(2,2);...
			 KI_BHA16(1,1)\KP_BHA16(1,1);  	KI_BHA16(2,2)\KP_BHA16(2,2);...
			 KI_CS03(1,1)\KP_CS03(1,1);  	KI_CS03(2,2)\KP_CS03(2,2);...
			 KI_BLT(1,1)\KP_BLT(1,1);  		KI_BLT(2,2)\KP_BLT(2,2);...
             KI_ITAE(1,1)\KP_ITAE(1,1);  	KI_ITAE(2,2)\KP_ITAE(2,2);...
             KI_GAR21(1,1)\KP_GAR21(1,1);  	KI_GAR21(2,2)\KP_GAR21(2,2)],'%.2f'),...
    num2str([0;  							0;...
             0;  							0;...
             KP_TW7(1,1)\KD_TW7(1,1);  	    KP_TW7(2,2)\KD_TW7(2,2);...
             KP_TW12(1,1)\KD_TW12(1,1);  	KP_TW12(2,2)\KD_TW12(2,2);...
			 KP_BHA16(1,1)\KD_BHA16(1,1);  	KP_BHA16(2,2)\KD_BHA16(2,2);...
			 0;  							0;...
			 0;  							0;...
             KP_ITAE(1,1)\KD_ITAE(1,1);  	KP_ITAE(2,2)\KD_ITAE(2,2);...
             KP_GAR21(1,1)\KD_GAR21(1,1);  	KP_GAR21(2,2)\KD_GAR21(2,2)],'%.2f'),...			 
	num2str([info_TW1(1,1).ITAE;			info_TW1(2,1).ITAE;...
             info_TW5(1,1).ITAE;			info_TW5(2,1).ITAE;...
             info_TW7(1,1).ITAE;			info_TW7(2,1).ITAE;...
             info_TW12(1,1).ITAE;			info_TW12(2,1).ITAE;...
			 info_BHA16(1,1).ITAE;			info_BHA16(2,1).ITAE;...
			 info_CS03(1,1).ITAE;			info_CS03(2,1).ITAE;...
             info_BLT(1,1).ITAE;			info_BLT(2,1).ITAE;...
             info_ITAE(1,1).ITAE;			info_ITAE(2,1).ITAE;...
             info_GAR21(1,1).ITAE;			info_GAR21(2,1).ITAE],'%.1f'),...
	num2str([info_TW1(1,2).ITAE;			info_TW1(2,2).ITAE;...
             info_TW5(1,2).ITAE;			info_TW5(2,2).ITAE;...
             info_TW7(1,2).ITAE;			info_TW7(2,2).ITAE;...
             info_TW12(1,2).ITAE;			info_TW12(2,2).ITAE;...
			 info_BHA16(1,2).ITAE;			info_BHA16(2,2).ITAE;...
			 info_CS03(1,2).ITAE;			info_CS03(2,2).ITAE;...
             info_BLT(1,2).ITAE;			info_BLT(2,2).ITAE;...
             info_ITAE(1,2).ITAE;			info_ITAE(2,2).ITAE;...
             info_GAR21(1,2).ITAE;			info_GAR21(2,2).ITAE],'%.1f'),...
	num2str([info_TW1(1,3).ITAE;			info_TW1(2,3).ITAE;...
             info_TW5(1,3).ITAE;			info_TW5(2,3).ITAE;...
             info_TW7(1,3).ITAE;			info_TW7(2,3).ITAE;...
             info_TW12(1,3).ITAE;			info_TW12(2,3).ITAE;...
			 info_BHA16(1,3).ITAE;			info_BHA16(2,3).ITAE;...
			 info_CS03(1,3).ITAE;			info_CS03(2,3).ITAE;...
             info_BLT(1,3).ITAE;			info_BLT(2,3).ITAE;...
             info_ITAE(1,3).ITAE;			info_ITAE(2,3).ITAE;...
             info_GAR21(1,3).ITAE;			info_GAR21(2,3).ITAE],'%.1f'),...
	num2str([info_TW1(1,1).TransientTime;	info_TW1(2,1).TransientTime;...
             info_TW5(1,1).TransientTime;	info_TW5(2,1).TransientTime;...
             info_TW7(1,1).TransientTime;	info_TW7(2,1).TransientTime;...
             info_TW12(1,1).TransientTime;	info_TW12(2,1).TransientTime;...
			 info_BHA16(1,1).TransientTime;	info_BHA16(2,1).TransientTime;...
			 info_CS03(1,1).TransientTime;	info_CS03(2,1).TransientTime;...
             info_BLT(1,1).TransientTime;	info_BLT(2,1).TransientTime;...
             info_ITAE(1,1).TransientTime;	info_ITAE(2,1).TransientTime;...
             info_GAR21(1,1).TransientTime;	info_GAR21(2,1).TransientTime],'%.1f'),...
	num2str([info_TW1(1,2).TransientTime;	info_TW1(2,2).TransientTime;...
             info_TW5(1,2).TransientTime;	info_TW5(2,2).TransientTime;...
             info_TW7(1,2).TransientTime;	info_TW7(2,2).TransientTime;...
             info_TW12(1,2).TransientTime;	info_TW12(2,2).TransientTime;...
			 info_BHA16(1,2).TransientTime;	info_BHA16(2,2).TransientTime;...
			 info_CS03(1,2).TransientTime;	info_CS03(2,2).TransientTime;...
             info_BLT(1,2).TransientTime;	info_BLT(2,2).TransientTime;...
             info_ITAE(1,2).TransientTime;	info_ITAE(2,2).TransientTime;...
             info_GAR21(1,2).TransientTime;	info_GAR21(2,2).TransientTime],'%.1f'),...
	num2str([info_TW1(1,3).TransientTime;	info_TW1(2,3).TransientTime;...
             info_TW5(1,3).TransientTime;	info_TW5(2,3).TransientTime;...
             info_TW7(1,3).TransientTime;	info_TW7(2,3).TransientTime;...
             info_TW12(1,3).TransientTime;	info_TW12(2,3).TransientTime;...
			 info_BHA16(1,3).TransientTime;	info_BHA16(2,3).TransientTime;...
			 info_CS03(1,3).TransientTime;	info_CS03(2,3).TransientTime;...
             info_BLT(1,3).TransientTime;	info_BLT(2,3).TransientTime;...
             info_ITAE(1,3).TransientTime;	info_ITAE(2,3).TransientTime;...
             info_GAR21(1,3).TransientTime;	info_GAR21(2,3).TransientTime],'%.1f')...
    );
tab5.Properties.RowNames={...
    'This Work #1 r-w_D',...
    'This Work #1 v-w_B',...
    'This Work #5 r-w_D',...
    'This Work #5 v-w_B',...    
    'This Work #7 r-w_D',...
    'This Work #7 v-w_B',...        
    'This Work #12 r-w_D',...
    'This Work #12 v-w_B',...            
    'BHA16 r-w_D',...
    'BHA16 v-w_B',...
    'CS03 r-w_D',...
    'CS03 v-w_B',...
    'BLT r-w_D',...
    'BLT v-w_B',...
    'ITAE r-w_D',...
    'ITAE v-w_B',...
    'GAR21 r-w_D',...
    'GAR21 v-w_B'};
tab5.Properties.VariableNames={'K_P','T_I','T_D','ITAE w_{D,r}','ITAE w_{B,r}','ITAE f','Sett T w_{D,r}','Sett T w_{B,r}','Sett T f'};
tab5
save('cases\WoodBerry\yalmip\results\data\yalmip_exWB_tb_rlt.mat','tab5','-append');
%--------------------------------------------------------------------------