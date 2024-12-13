clear, close all
%--------------------------------------------------------------------
%SANTANDER et al. (2022) FCC MODEL ----------------------------------
% xfcc=[1564.06758445685,616,024.0440758965397,24.9000000145994,40.0420669348231,28.0000000041903,3.16075658859897,35.0875282505755,0.0101439627122109,0.00250580481493606,1249.99997383664,9350.79072581430,271210.435908199,246.810877809621,98098.9991960623,968.999859177516,2.68045714195376,30.4463683645626,14.6380093682365,-24.2765676899209,-0.612986761842691,0.611929610398922,-5641.47052429781,-30153150.1669079,-1879.52204006479,0.00468987379029097,0.0481251434795964,0.527273962836787,0.357213788281919,2.47339657486076,1.41455228053760,2.70179775410060,2.44381704476869,0.631637205219422,0.627951527686494,457.718461235902,6.31203102548586];
% ufcc=[165,eps,1,eps,eps,eps,0.431260000000000,eps,eps,616,24.9000000000000,28,1250,98100,969];% LAST 6 COMPONENTS ARE THE SET POINTS FOR FCC. 
% dist=[75,25,460.900000000000,0.900000000000000]; %DISTURBANCES
% Flpg=457.718467816417;
% Tcondenser=310.707095521277;
% 
% xfcc=xfcc(:); ufcc=ufcc(:); dist=dist(:);
% [~,yp]=FCC(xfcc,dist,ufcc,Flpg,Tcondenser,1);
% yp=yp(:);

load('cases\FCC\nominal\dynamic_ss.mat');
xfcc=xfcc(:);
ufcc=ufcc(:); ufcc(find(ufcc==0))=eps;
dist=dist(:);
[~,yp]=FCC(xfcc,dist,ufcc,Flpg,Tcondenser,1);
yp=yp(:);

x0=xfcc([(1:2),(5:19),(26:35),37]);
u0=[ufcc(1:9);[yp(23)*100/118.75;yp(21);yp(20)*160/100;yp(22);yp(24)];dist(1:3)];

sys=sfunc2ss(@FCC_sfunc,x0,u0,[]);

idu=(10:14)';
idw=[(1:9)';(15:17)'];

idu0=unique([idu;idw]);
for i=1:size(idu,1); idut(i)=find(idu(i)==idu0); end; idut=idut(:);
for i=1:size(idw,1); idwt(i)=find(idw(i)==idu0); end; idwt=idwt(:);

%T2,  P6,  Treg, Wr,   Tr
%(5), (3), (7),  (28), (6) 
idy=[3;5;6;7;28];
idy0=idy; idyt=idy0;

%Control action structure matrix
%	 P6 T2	Tr	Treg	Wr
%idy:3	5	6	7		28  idu:
Tk=	[0	1	0	0		0	%10 p1a (Ff IN MANUSCRIPT)
	 -1	0	0	0		0	%11 V7a (Ffg IN MANUSCRIPT)
	 0	0	0	1		0	%12 p6a (Fa IN MANUSCRIPT)
	 0	0	0	0		-1	%13 V3a (Fsc IN MANUSCRIPT)
	 0	0	1	0		0];	%14 V2a (Frgc IN MANUSCRIPT)

%Instrumentation gains
Kmy=diag([...
	1e3,	...y(3):  	P6		Regenerator pressure
	1,      ...y(5):  	T2   	Preheater outlet temperature
	1,      ...y(6):  	Tr  	Reactor riser temperature
	1,  	...y(7): 	Treg    Regenerator temperature
	1e-4]   ...y(28): 	Wr      Inventory catalyst in reactor
	);

%--------------------------------------------------------------------
%OPEN-LOOP SYSTEM MATRICES FOR CL VARIABLES WITH PI CONTROL ---------	 
A=sys.A;
B=sys.B(:,idu0);
Cs=sys.C(idyt,:);
Cy=Kmy*Cs;
Ds=sys.D;

%--------------------------------------------------------------------
%MODEL VALIDATION ---------------------------------------------------
%check system eigenvalues
deig=eig(A);

%Need to check if Dyu and Dyw are null matrices
Bu=B(:,idut);
Bw=B(:,idwt);
Dyu=Kmy*Ds(idyt,idut);
Dyw=Kmy*Ds(idyt,idwt);
%--------------------------------------------------------------------
%AUGMENTED PI OPEN-LOOP SYSTEM MATRICES -----------------------------
Nx=size(A,1);   %State variables size
Ny=size(Cy,1);  %Controlled output variables size
Nu=size(Bu,2);  %Input variables size
Nw=size(Bw,2);  %Exogenous input variables size

%Open-loop autonomous system matrix (Nx+Ny,Nx+Ny)
A_=	[A              zeros(Nx,Ny)
	 Cy           	zeros(Ny)];
 
%Open-loop system input gain matrix (Nx+Ny,Nu)
Bu_=[Bu
	 zeros(Ny,Nu)];
 
%Open-loop system exogenous input gain matrix (Nx+Ny,Nw)
Bw_=[Bw
	 zeros(Ny,Nw)];		 

%Output-State matrix relationship (Ny+Ny,Nx+Ny)
Cs_=[Cs          	zeros(Ny)
	 zeros(Ny,Nx)   eye(Ny)];
 
%Output-State matrix relationship (Ny+Ny,Nx+Ny)
Cy_=[Cy          		zeros(Ny)
	zeros(Ny,Nx)   	eye(Ny)];
%--------------------------------------------------------------------
%POLYTOPIC LDI SYSTEM -----------------------------------------------
Poly.Nl=1; Poly.Nx=Nx; Poly.Ny=Ny; Poly.Nu=Nu; Poly.Nw=Nw;
Poly.A{1}=A;   Poly.A0=A;
Poly.Bu{1}=Bu; Poly.Bu0=Bu;
Poly.Bw{1}=Bw; Poly.Bw0=Bw;
Poly.Cy{1}=Cy; Poly.Cy0=Cy;
%--------------------------------------------------------------------
%LMI PROBLEM TATEMENT ----------------------------------------------
%LMI problem parameters

%Proportional gain nominal values
KPnom=diag([...
	1.5;...1500 ...y(3):  	P6		Regenerator pressure
	5; 			...y(5):  	T2   	Preheater outlet temperature
	0.5; 		...y(6):  	Tr  	Reactor riser temperature
	1.5; 		...y(7): 	Treg    Regenerator temperature
	1.0]...1e-4 ...y(28): 	Wr     	Inventory catalyst in reactor
	);

%Proportional gain upper and lower bounds
dlt=0.5;
KPlb=dlt*KPnom;
KPub=dlt\KPnom;
	
%Integral time nominal values
TInom=diag([...
	200; 		...y(3):  	P6		Regenerator pressure
	100; 		...y(5):  	T2   	Preheater outlet temperature
	200; 		...y(6):  	Tr  	Reactor riser temperature
	800; 		...y(7): 	Treg    Regenerator temperature
	1000] 		...y(28): 	Wr      Inventory catalyst in reactor
	);
	
%Integral time upper and lower bounds
dlt=0.5;
TIlb=dlt*TInom;
TIub=dlt\TInom;

% Use when not tuning by LMI SOF stabilization
% KPpi=(KPnom*Tk); KIpi=TInom\(KPnom*Tk);

myfile='cases\FCC\lmilab_FCC_SOFS_YALMIP_PI';
diary(myfile);
%--------------------------------------------------------------------
%PI RESULTS ---------------------------------------------------------
%PI tuning by SOF stabilization
timerVal = tic;
[KPpi,KIpi,~,TIpi,~,alppi,P_pi,Sgm_FPpi,epsl_FPpi,~,~]=...
	PID_SOFS(Poly,Tk,[],KPlb,KPub,TIlb,TIub,[],[],'PI');
elapsedTime_pi = toc(timerVal);

timerVal = tic;
[KPpi_ILMI,KIpi_ILMI,~,TIpi_ILMI,~,alppi_ILMI,P_pi_ILMI,iterpi_ILMI]=...
    PID_SOFS_ILMI(Poly,Tk,[],KPlb,KPub,TIlb,TIub,[],[],'PI');
elapsedTime_pi_ILMI = toc(timerVal);
%--------------------------------------------------------------------
diary off

%Output step response for PI control
tstart=0; tend=60*500/2; nd=1000; t=(tstart:(tend-tstart)/nd:tend)'; %time domain
Nt=size(t,1);

%Closed-loop system input gain matrix (Nx+Ny,Nu)
Br_=[Bu*KPpi*Kmy
	 -eye(Ny,Nu)*Kmy];
Br_ILMI=[Bu*KPpi_ILMI*Kmy
	 -eye(Ny,Nu)*Kmy];

%Santander et al. (2022) set-point tracking simulation
sp=[2.5;    ...y(5):  	T2   	Preheater outlet temperature
	2];     ...y(6):  	Tr  	Reactor riser temperature

y0=Cs*x0;

syscl_SP=ss(A_-Bu_*[KPpi,KIpi]*Cy_, Br_(:,[2;3]), Cs_(1:Ny,:), []);
syscl_SP_ILMI=ss(A_-Bu_*[KPpi_ILMI,KIpi_ILMI]*Cy_, Br_ILMI(:,[2;3]), Cs_(1:Ny,:), []);
stepOpt=stepDataOptions('StepAmplitude',sp);
figure(1); step(syscl_SP,t,stepOpt); hold on; step(syscl_SP_ILMI,t,stepOpt)
ycl_SP=step(syscl_SP,t,stepOpt); ycl_SP_ILMI=step(syscl_SP_ILMI,t,stepOpt);
ycl_SP=ones(Nt,1)*y0'+ycl_SP; ycl_SP_ILMI=ones(Nt,1)*y0'+ycl_SP_ILMI; 
syscl_SP_Hinf=hinfnorm(syscl_SP);
syscl_SP_ILMI_Hinf=hinfnorm(syscl_SP_ILMI);
%--------------------------------------------------------------------
%Output step response for PI control
%Santander et al. (2022) disturbance rejection simulation
w=-5;...  Feed quality disturbance

syscl_DR=ss(A_-Bu_*[KPpi,KIpi]*Cy_, Bw_(:,find(idw==16)), Cs_(1:Ny,:), []);
syscl_DR_ILMI=ss(A_-Bu_*[KPpi_ILMI,KIpi_ILMI]*Cy_, Bw_(:,find(idw==16)), Cs_(1:Ny,:), []);
stepOpt=stepDataOptions('StepAmplitude',w);
figure(2); step(syscl_DR,t,stepOpt); hold on; step(syscl_DR_ILMI,t,stepOpt)
ycl_DR=step(syscl_DR,t,stepOpt); ycl_DR_ILMI=step(syscl_DR_ILMI,t,stepOpt);
ycl_DR=ones(Nt,1)*y0'+ycl_DR; ycl_DR_ILMI=ones(Nt,1)*y0'+ycl_DR_ILMI;
syscl_DR_Hinf=hinfnorm(syscl_DR);
syscl_DR_H_ILMI_Hinf=hinfnorm(syscl_DR_ILMI);
%--------------------------------------------------------------------
save(strjoin({myfile,'.mat'},''));
