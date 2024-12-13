clear
%--------------------------------------------------------------------
%WOOD & BERRY (1973) MODEL ------------------------------------------
%Transfer function model
s=tf('s');
dl=[[1 3 8.1];[7 3 3.4]];
G(1,1)=12.8*exp(-dl(1,1)*s)/(16.7*s+1);     %wD(s)/r(s)
G(1,2)=-18.9*exp(-dl(1,2)*s)/(21*s+1);      %wD(s)/v(s)
G(2,1)=6.6*exp(-dl(2,1)*s)/(10.9*s+1);      %wB(s)/r(s)
G(2,2)=-19.4*exp(-dl(2,2)*s)/(14.4*s+1);    %wB(s)/v(s)
G(1,3)=3.8*exp(-dl(1,3)*s)/(14.9*s+1);      %wD(s)/f(s)
G(2,3)=4.9*exp(-dl(2,3)*s)/(13.2*s+1);      %wB(s)/f(s)
num=get(G,'num');
den=get(G,'den');

%Wood & Berry (1973) ss model neglecting input delays
sys=ss;
sys.A=-inv(diag([den{1,1}(1);den{1,2}(1);den{1,3}(1);den{2,1}(1);den{2,2}(1);den{2,3}(1)]));
sys.B=(-sys.A)*...
    [num{1,1}(2) 0              0
     0           num{1,2}(2)    0
     0           0              num{1,3}(2)
     num{2,1}(2) 0              0
     0           num{2,2}(2)    0
     0           0              num{2,3}(2)];
sys.C=[[1 0];[1 0];[1 0];[0 1];[0 1];[0 1]]';
sys.D=[];

Nx=size(sys.A,1); Ny=size(sys.C,1); Nu=size(sys.B,2);
Nsd=5;	%Number of discretizations

%Reflux 'r' input delay virtual ss model
Trs=max(dl(:,1))/Nsd;
sysr=ss;
sysr.A=(diag(-1*ones(1,Nsd))+diag(ones(1,Nsd-1),-1))/Trs;
sysr.B=zeros(Nsd,Nu); sysr.B(1,1)=1/Trs;
sysr.C=zeros(Ny,Nsd);
sysr.D=[];

%Bottom vapor 'v' input delay virtual ss model
Tvs=max(dl(:,2))/Nsd;
sysv=ss;
sysv.A=(diag(-1*ones(1,Nsd))+diag(ones(1,Nsd-1),-1))/Tvs;
sysv.B=zeros(Nsd,Nu); sysv.B(1,2)=1/Tvs;
sysv.C=zeros(Ny,Nsd);
sysv.D=[];

%feed 'f' input delay virtual ss model
Tfs=max(dl(:,3))/Nsd;
sysf=ss;
sysf.A=(diag(-1*ones(1,Nsd))+diag(ones(1,Nsd-1),-1))/Tfs;
sysf.B=zeros(Nsd,Nu); sysf.B(1,3)=1/Tfs;
sysf.C=zeros(Ny,Nsd);
sysf.D=[];

%Delay time Continuos-time-approximation (CTA) ss model
sysm=ss;

%Autonomous part
sysm.A=zeros(Nx+Nsd+Nsd+Nsd);
sysm.A(1:Nsd,1:Nsd)=sysr.A;
sysm.A(Nsd+1:Nsd+Nsd,Nsd+1:Nsd+Nsd)=sysv.A;
sysm.A(Nsd+Nsd+1:Nsd+Nsd+Nsd,Nsd+Nsd+1:Nsd+Nsd+Nsd)=sysf.A;
sysm.A(Nsd+Nsd+Nsd+1:end,Nsd+Nsd+Nsd+1:end)=sys.A;

sysm.B=[sysr.B;sysv.B;sysf.B;zeros(6,Nu)];  %State-Input
sysm.C=[sysr.C,sysv.C,sysf.C,sys.C];        %Output-State
sysm.D=[];                                  %Output-Input

%Feeding to real states the delayed inputs as virtual states
sysm.A(Nsd+Nsd+Nsd+1,round(dl(1,1)/Trs))=sys.B(1,1);
sysm.A(Nsd+Nsd+Nsd+2,Nsd+round(dl(1,2)/Tvs))=sys.B(2,2);
sysm.A(Nsd+Nsd+Nsd+3,Nsd+Nsd+round(dl(1,3)/Tfs))=sys.B(3,3);
sysm.A(Nsd+Nsd+Nsd+4,round(dl(2,1)/Trs))=sys.B(4,1);
sysm.A(Nsd+Nsd+Nsd+5,Nsd+round(dl(2,2)/Tvs))=sys.B(5,2);
sysm.A(Nsd+Nsd+Nsd+6,Nsd+Nsd+round(dl(2,3)/Tfs))=sys.B(6,3);

%Uncomment the line below for LTI CTA model validation
%figure(1); step(G); hold on; step(sysm); close figure(1);
%--------------------------------------------------------------------
%OPEN-LOOP SYSTEM MATRICES ------------------------------------------
A=sysm.A;
Bu=sysm.B(:,1:Ny);
Bw=sysm.B(:,Ny+1:end);
Cy=sysm.C;

Nx=size(A,1);   %State variables size
Ny=size(Cy,1);  %Controlled output variables size
Nu=size(Bu,2);  %Input variables size
Nw=size(Bw,2);  %Exogenous input variables size

PhiD=6e1*eye(Ny); %Derivative kick filter
% PhiD=6e2*eye(Ny); %Derivative kick filter
%--------------------------------------------------------------------
%{
%POLYTOPIC LDI SYSTEM -----------------------------------------------
%Nominal LTI system
Poly.A0=A;
Poly.Bu0=Bu;
Poly.Bw0=Bw;
Poly.Cy0=Cy;

dl0=dl; Nl=0; dl={};
for i=1:size(dl0,1)
    for j=1:size(dl0,2)
        dl{Nl+1}=zeros(size(dl0,1),size(dl0,2));
        dl{Nl+1}(i,j)=dl0(i,j)*0.1;
        Nl=Nl+1;
    end
end
for i=1:Nl
	dl{i}=dl0+dl{i};
end
dl=[dl0,dl]; Nl=Nl+1;
Poly.Nl=Nl; Poly.Nx=Nx; Poly.Ny=Ny; Poly.Nw=Nw;
Poly.A=cell(Nl,1);
Poly.Bu=cell(Nl,1);
Poly.Bw=cell(Nl,1);
Poly.Cy=cell(Nl,1);

for i=1:Nl
Trs=max(dl{i}(:,1))/Nsd; Tvs=max(dl{i}(:,2))/Nsd; Tfs=max(dl{i}(:,3))/Nsd;

%Autonomous part
Poly.A{i}=zeros(size(sys.A,1)+Nsd+Nsd+Nsd);
Poly.A{i}(1:Nsd,1:Nsd)=(diag(-1*ones(1,Nsd))+diag(ones(1,Nsd-1),-1))/Trs;
Poly.A{i}(Nsd+1:Nsd+Nsd,Nsd+1:Nsd+Nsd)=(diag(-1*ones(1,Nsd))+diag(ones(1,Nsd-1),-1))/Tvs;
Poly.A{i}(Nsd+Nsd+1:Nsd+Nsd+Nsd,Nsd+Nsd+1:Nsd+Nsd+Nsd)=(diag(-1*ones(1,Nsd))+diag(ones(1,Nsd-1),-1))/Tfs;
Poly.A{i}(Nsd+Nsd+Nsd+1:end,Nsd+Nsd+Nsd+1:end)=sys.A;
Poly.A{i}(Nsd+Nsd+Nsd+1,round(dl{i}(1,1)/Trs))=sys.B(1,1);
Poly.A{i}(Nsd+Nsd+Nsd+2,Nsd+round(dl{i}(1,2)/Tvs))=sys.B(2,2);
Poly.A{i}(Nsd+Nsd+Nsd+3,Nsd+Nsd+round(dl{i}(1,3)/Tfs))=sys.B(3,3);
Poly.A{i}(Nsd+Nsd+Nsd+4,round(dl{i}(2,1)/Trs))=sys.B(4,1);
Poly.A{i}(Nsd+Nsd+Nsd+5,Nsd+round(dl{i}(2,2)/Tvs))=sys.B(5,2);
Poly.A{i}(Nsd+Nsd+Nsd+6,Nsd+Nsd+round(dl{i}(2,3)/Tfs))=sys.B(6,3);

%State-Input
rB=[[1/Trs,zeros(1,Ny-1)];zeros(Nsd-1,Ny)];
vB=[[zeros(1,Ny-1),1/Tvs];zeros(Nsd-1,Ny)];
Poly.Bu{i}=[rB;vB;zeros(Nsd,Ny);zeros(6,Ny)];

%State-Exogenous input
fB=[1/Tfs;zeros(Nsd-1,Nw)];
Poly.Bw{i}=[zeros(Nsd,Nw);zeros(Nsd,Nw);fB;zeros(6,Nw)];	

%Output-State
Poly.Cy{i}=Cy;

%System eigenvalues
Poly.deig{i}=eig(Poly.A{i});
end

%{
%Uncommenting this section for PLDI CTA model validation
for i=1:Nl
G(1,1)=12.8*exp(-dl{i}(1,1)*s)/(16.7*s+1);     %wD(s)/r(s)
G(1,2)=-18.9*exp(-dl{i}(1,2)*s)/(21*s+1);      %wD(s)/v(s)
G(2,1)=6.6*exp(-dl{i}(2,1)*s)/(10.9*s+1);      %wB(s)/r(s)
G(2,2)=-19.4*exp(-dl{i}(2,2)*s)/(14.4*s+1);    %wB(s)/v(s)
G(1,3)=3.8*exp(-dl{i}(1,3)*s)/(14.9*s+1);      %wD(s)/f(s)
G(2,3)=4.9*exp(-dl{i}(2,3)*s)/(13.2*s+1);      %wB(s)/f(s)
figure(i); step(G); hold on; step(ss(Poly.A{i},[Poly.Bu{i},Poly.Bw{i}],Poly.Cy{i},[]));
end
%close all
%}
%}

%Uncommenting this section for LTI system instead of PLDI
clear Poly
Poly.A0=A;
Poly.Bu0=Bu;
Poly.Bw0=Bw;
Poly.Cy0=Cy;
Poly.Nl=1; Poly.Nx=Nx; Poly.Ny=Ny; Poly.Nu=Nu; Poly.Nw=Nw;
Poly.A{1}=Poly.A0;
Poly.Bu{1}=Poly.Bu0;
Poly.Bw{1}=Poly.Bw0;
Poly.Cy{1}=Poly.Cy0;
Poly.deig{1}=eig(Poly.A0);
%}
%-------------------------------------------------------------------- 
%LMI PROBLEM STATEMENT ----------------------------------------------
%control action structure matrix
Tk=...
	[1  	0
     0     -1];

%Control parameters sector constratins
%Constraints with sluggish result
%{
KPlb=...                %Proportional gain lower bound
    [0.1   0
     0    -0.1];
KPub=...                %Proportional gain upper bound
    [0.5    0
     0    	-0.5];
TIlb=...                %Integral time lower bound
    [5  	0
     0    	5]; 
TIub=...                %Integral time upper bound
    [300  	0
     0    	300];
TDlb=...                %Derivative time lower bound
    [0  	0
     0    	0];
TDub=...                %Derivative time upper bound
    [1      0
     0    	1];	 
%}

%Better constraints after attempt with below
KPlb=...                %Proportional gain lower bound
    [1e-4   0
     0      -1e-4];
KPub=...                %Proportional gain upper bound
    [0.1    0
     0    	-0.1];
TIlb=...                %Integral time lower bound
    [5  	0
     0    	5];
TIub=...                %Integral time upper bound
    [50  	0
     0    	50];
TDlb=...                %Derivative time lower bound
    [1e-3  	0
     0    	1e-3];
TDub=...                %Derivative time upper bound
    [1      0
     0    	1];	 

%Diagonalize KP, TI and TD upper and lower bounds
KPlb=Tk*KPlb; KPub=Tk*KPub;
TIlb=abs(Tk)*TIlb; TIub=abs(Tk)*TIub;
TDlb=abs(Tk)*TDlb; TDub=abs(Tk)*TDub;

%Upper bound on L2 norm
gL2=1e3;
% gL2=1e2;
% gL2=1e1;

% myfile=strjoin({'cases\WoodBerry\yalmip\yalmip_exWB_SOFS_L2N_NOTITER_Nd=',num2str(Nsd),'_Nl=',num2str(Poly.Nl),'_PhiD=',num2str(nthroot(det(PhiD),Ny)),'_gL2=',num2str(gL2)},'');
myfile=strjoin({'cases\WoodBerry\yalmip\yalmip_exWB_SOFS_L2N_ILMI_Nd=',num2str(Nsd),'_Nl=',num2str(Poly.Nl),'_PhiD=',num2str(nthroot(det(PhiD),Ny)),'_gL2=',num2str(gL2)},'');
diary(myfile);
%--------------------------------------------------------------------
%AUGMENTED PI OPEN-LOOP SYSTEM MATRICES -----------------------------
%Open-loop autonomous system matrix (Nx+Ny,Nx+Ny)
A_pi= [A              zeros(Nx,Ny)
     Cy             zeros(Ny,Ny)];

%Open-loop system input gain matrix (Nx+Ny,Nu)
Bu_pi=[Bu
	 zeros(Ny,Nu)];
	 
%Open-loop system exogenous input gain matrix (Nx+Ny,Nw)
Bw_pi=[Bw
	 zeros(Ny,Nw)];

%Controlled Output-State matrix relationship (Ny+Ny,Ny+Ny)
Cy_pi=[Cy             zeros(Ny,Ny)
	 zeros(Ny,Nx)   eye(Ny,Ny)];
%--------------------------------------------------------------------
%PI RESULTS ---------------------------------------------------------
%PI tuning by SOF stabilization
disp('----- start of PI calculation ---------------------------------------');
% timerVal = tic;
% [KPpi,KIpi,~,TIpi,~,alppi,P_pi,Sgm_FPpi,epsl_FPpi,~,~]=...
%     PID_SOFS_L2N_YALMIP(Poly,gL2,Tk,[],KPlb,KPub,TIlb,TIub,[],[],'PI');
% elapsedTime_pi = toc(timerVal);

timerVal = tic;
[KPpi_ILMI,KIpi_ILMI,~,TIpi_ILMI,~,alppi_ILMI,P_pi_ILMI,iterpi_ILMI]=...
    PID_SOFS_L2N_YALMIP_ILMI(Poly,gL2,Tk,[],KPlb,KPub,TIlb,TIub,[],[],'PI');
elapsedTime_pi_ILMI = toc(timerVal);
disp('----- end of PI calculation -----------------------------------------');
%--------------------------------------------------------------------
%AUGMENTED PID OPEN-LOOP SYSTEM MATRICES ----------------------------
%Open-loop autonomous system matrix (Nx+2*Ny,Nx+2*Ny)
A_pid= [A           zeros(Nx,Ny)	zeros(Nx,Ny)
		Cy          zeros(Ny)		zeros(Ny)
		PhiD*Cy*A 	zeros(Ny)		-PhiD*eye(Ny)];

%Open-loop system input gain matrix (Nx+2*Ny,Nu)
Bu_pid=[Bu
		zeros(Ny,Nu)
		PhiD*Cy*Bu];
	 
%Open-loop system exogenous input gain matrix (Nx+2*Ny,Nw)
Bw_pid=[Bw
		zeros(Ny,Nw)
		PhiD*Cy*Bw];

%Controlled Output-State matrix relationship (Ny+2*Ny,Ny+2*Ny)
Cy_pid=[Cy             	zeros(Ny)       zeros(Ny)
		zeros(Ny,Nx)   	eye(Ny)			zeros(Ny)
		zeros(Ny,Nx)	zeros(Ny)		eye(Ny)];
	 
%-------------------------------------------------------------------- 
%PID RESULTS ---------------------------------------------------------
% PID tuning by SOF stabilization
disp('----- start of PID calculation --------------------------------------');
% timerVal = tic;
% [KPpid,KIpid,KDpid,TIpid,TDpid,alppid,P_pid,Sgm_FPpid,epsl_FPpid,Sgm_Pprppid,epsl_Pprppid]=...
%     PID_SOFS_L2N_YALMIP(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,'PID');	
% elapsedTime_pid = toc(timerVal);

timerVal = tic;
[KPpid_ILMI,KIpid_ILMI,KDpid_ILMI,TIpid_ILMI,TDpid_ILMI,alppid_ILMI,P_pid_ILMI,iterpid_ILMI]=...
    PID_SOFS_L2N_YALMIP_ILMI(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,'PID');
elapsedTime_pid_ILMI = toc(timerVal);
disp('----- end of PID calculation ----------------------------------------');
%--------------------------------------------------------------------
diary off;
save(myfile);

%Upper bound on L2 norm
% gL2=1e3;
gL2=1e2;
% gL2=1e1;

% myfile=strjoin({'cases\WoodBerry\yalmip\yalmip_exWB_SOFS_L2N_NOTITER_Nd=',num2str(Nsd),'_Nl=',num2str(Poly.Nl),'_PhiD=',num2str(nthroot(det(PhiD),Ny)),'_gL2=',num2str(gL2)},'');
myfile=strjoin({'cases\WoodBerry\yalmip\yalmip_exWB_SOFS_L2N_ILMI_Nd=',num2str(Nsd),'_Nl=',num2str(Poly.Nl),'_PhiD=',num2str(nthroot(det(PhiD),Ny)),'_gL2=',num2str(gL2)},'');
diary(myfile);
%PI RESULTS ---------------------------------------------------------
%PI tuning by SOF stabilization
disp('----- start of PI calculation ---------------------------------------');
% timerVal = tic;
% [KPpi,KIpi,~,TIpi,~,alppi,P_pi,Sgm_FPpi,epsl_FPpi,~,~]=...
%     PID_SOFS_L2N_YALMIP(Poly,gL2,Tk,[],KPlb,KPub,TIlb,TIub,[],[],'PI');
% elapsedTime_pi = toc(timerVal);

timerVal = tic;
[KPpi_ILMI,KIpi_ILMI,~,TIpi_ILMI,~,alppi_ILMI,P_pi_ILMI,iterpi_ILMI]=...
    PID_SOFS_L2N_YALMIP_ILMI(Poly,gL2,Tk,[],KPlb,KPub,TIlb,TIub,[],[],'PI');
elapsedTime_pi_ILMI = toc(timerVal);
disp('----- end of PI calculation -----------------------------------------');
%--------------------------------------------------------------------
%PID RESULTS ---------------------------------------------------------
% PID tuning by SOF stabilization
disp('----- start of PID calculation --------------------------------------');
% timerVal = tic;
% [KPpid,KIpid,KDpid,TIpid,TDpid,alppid,P_pid,Sgm_FPpid,epsl_FPpid,Sgm_Pprppid,epsl_Pprppid]=...
%     PID_SOFS_L2N_YALMIP(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,'PID');	
% elapsedTime_pid = toc(timerVal);

timerVal = tic;
[KPpid_ILMI,KIpid_ILMI,KDpid_ILMI,TIpid_ILMI,TDpid_ILMI,alppid_ILMI,P_pid_ILMI,iterpid_ILMI]=...
    PID_SOFS_L2N_YALMIP_ILMI(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,'PID');
elapsedTime_pid_ILMI = toc(timerVal);
disp('----- end of PID calculation ----------------------------------------');
%--------------------------------------------------------------------
diary off;
save(myfile);

%Upper bound on L2 norm
% gL2=1e3;
% gL2=1e2;
gL2=1e1;

% myfile=strjoin({'cases\WoodBerry\yalmip\yalmip_exWB_SOFS_L2N_NOTITER_Nd=',num2str(Nsd),'_Nl=',num2str(Poly.Nl),'_PhiD=',num2str(nthroot(det(PhiD),Ny)),'_gL2=',num2str(gL2)},'');
myfile=strjoin({'cases\WoodBerry\yalmip\yalmip_exWB_SOFS_L2N_ILMI_Nd=',num2str(Nsd),'_Nl=',num2str(Poly.Nl),'_PhiD=',num2str(nthroot(det(PhiD),Ny)),'_gL2=',num2str(gL2)},'');
diary(myfile);
%--------------------------------------------------------------------
%PI RESULTS ---------------------------------------------------------
%PI tuning by SOF stabilization
disp('----- start of PI calculation ---------------------------------------');
% timerVal = tic;
% [KPpi,KIpi,~,TIpi,~,alppi,P_pi,Sgm_FPpi,epsl_FPpi,~,~]=...
%     PID_SOFS_L2N_YALMIP(Poly,gL2,Tk,[],KPlb,KPub,TIlb,TIub,[],[],'PI');
% elapsedTime_pi = toc(timerVal);

timerVal = tic;
[KPpi_ILMI,KIpi_ILMI,~,TIpi_ILMI,~,alppi_ILMI,P_pi_ILMI,iterpi_ILMI]=...
    PID_SOFS_L2N_YALMIP_ILMI(Poly,gL2,Tk,[],KPlb,KPub,TIlb,TIub,[],[],'PI');
elapsedTime_pi_ILMI = toc(timerVal);
disp('----- end of PI calculation -----------------------------------------');
%--------------------------------------------------------------------
%PID RESULTS ---------------------------------------------------------
% PID tuning by SOF stabilization
disp('----- start of PID calculation --------------------------------------');
% timerVal = tic;
% [KPpid,KIpid,KDpid,TIpid,TDpid,alppid,P_pid,Sgm_FPpid,epsl_FPpid,Sgm_Pprppid,epsl_Pprppid]=...
%     PID_SOFS_L2N_YALMIP(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,'PID');	
% elapsedTime_pid = toc(timerVal);

timerVal = tic;
[KPpid_ILMI,KIpid_ILMI,KDpid_ILMI,TIpid_ILMI,TDpid_ILMI,alppid_ILMI,P_pid_ILMI,iterpid_ILMI]=...
    PID_SOFS_L2N_YALMIP_ILMI(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,'PID');
elapsedTime_pid_ILMI = toc(timerVal);
disp('----- end of PID calculation ----------------------------------------');
%--------------------------------------------------------------------
diary off;
save(myfile);
%}