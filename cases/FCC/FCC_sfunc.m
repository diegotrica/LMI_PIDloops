function [sys,x0,str,ts,simStateCompliance]=FCC_sfunc_v2(t,x,u,flag,cinit,par)
%FCC_sfunc - Fluid catalytic cracker unit model S-function
%
%   This model describes the dynamic behavior the fluid catalytic 
%	cracker (FCC) of Santander et al (2022).
%
%	[later type here a brief description of the adaptations made to 
%	build this S-function...]
%
%   [sys,x0,str,ts,simStateCompliance]=FCC_sfunc(t,x,u,flag,cinit,par)
%
%   sys is the S-function output, which depends on which flag is 
%   being inputted by simulink. x0 are the states initial condi-
%   tions. See the S-function Level-1 MATLAB reference page for more
%   details about str, ts and simStateCompliance S-Function ouputs.
%
%   t is the independent variable time. x are the state variables.
%   u are the input variables, combined into one vector of manipula-
%   ted and disturbs. flag is the Simulink variable inputted to get a
%   specific output from the S-function. par are the model parameters.
%
%   2024. Coded by Diego Trica (diegotrica@gmail.com)
%
%	References:
%	[1] Santander O, Kuppuraj V, Harrison CA, Baldea M. An open source
%	fluid catalytic cracker - fractionator model to support the 
%	development and benchmarking of process control, machine learning
%	and operation strategies, Computers and Chemical Engineering, 
%	2022, 164, DOI: 10.1016/j.compchemeng.2022.107900.

%(NEW) GLOBAL VARIABLES EXCHANGED BETWEEN SUBMODELS DURING SIMULATION
global Fair FH vs rhocdl Crgc zbed P6 epse Treg Fsg
global Tcyc o2cyc cocyc co2cyc
global xo2 xco xco2
global Wris Tr Frgc P4

%>>>>>>>>>Define distubances
	dist=u(15:end);
	if size(dist,1)>4
		error('dist must be a 4x1 vector');
	end
	Tatm = dist(1); %Ambient temperature disturbance
	API = dist(2); %Feed quality disturbance
	T1   = dist(3); %Feed temperature disturbance
	APInominal=25; %Nominal API (APInominal)
   
%'>>>>>>>Feed System Constants' 
   taufb=200.0; 
   DHfu=1000.0*(API/APInominal)^0.5; 
   UAf=25.0;
   a1=0.15;     
   a2=200.0;    
   taufo=60.0;
   F5nom=34; %Nominal Flow of fuel to furnace

%(NEW) COMPONENT MOLAR MASS FOR FRACTIONATOR SUBMODEL WHICH IS NOT CONSIDERED
%>>>>>>>>Component molar mass.'

    % MWp5=446.094358949534;%(g/mol)
    % MWp4=378.894679684485;
    % MWp3=292.185529267655;
    % MWp2=206.951432493898;
    % MWp1=120.941794467086;
    % MWc5=85.1347950005991;
    % MWb=58.12;
    % MWp=44.1;
    % MWe=30.07;
    % MWm=16.04;

%(NEW) CONSTANTS FOR WET GAS COMPRESSOR WHICH IS NOT CONSIDERED
 %>>>>>>>>Wet Gas Compressor Constants.'

	% k13=0.01;	
	% k11=1200.0212842935; 
	% Pvru=75.0;	

 %>>>>>>>>Lift Air BLower Constants.'

	vslip=2.2;	
	taufil=40.0;	
	hlift=34.0;
	sb=5950.0;	
	samin=5000.0;	
	kavg=1.39;
	etapla=1.0;	
	Tdla=225.0;		
	k9=10.0;
	Vdla=200.0;	
	klift=5.0;		
	k8=5.0;
    
 %>>>>>>>>Combustion Air Blower Constants.'

	Tdca=190.0;
	kca=40.0;		
	k6=250.0;	
	k7=15.0;
	Vcms=200.0;	
	Vdca=1000.0;

 %>>>>>>>>Regenerator Constants.'

	rgarea=590.0; 
	sparea=7.0;		
	cpCO2=11.0;
	cpair=7.08;		
	cpc=0.31;		
	cpCO=7.28;
	cpH2O=8.62;		
	cpN2=7.22;		
	cpO2=7.62; 
	k14=1.3464; 			
	MI=200000.0;	
	Qe=556.0;
	zcyc=45.0;		
	zsp=13.0;		
	R=10.73;					
	zlp=11.0;
	rhoprt=68.0;		
	rhoc=45.0;		
	Patm=14.7;		
	delH1=46368.0;
	fof=424;			
	hsp=20.0;		
	Tair=270.0;
	delH2=169080.0;	
	delhH=60960.0;	
	Tf=459.6;
    
 %>>>>>>Catalyst Recirculation:'

	areaur=3.7;  
	areaus=5.2;     
	Astrp=60.0;   
	Alp=8.73;    
	Lur=56.0;    
	Lus=56.0;       
	Etap=155.0;
	Eloi=124.5;  
	Estrp=130.0;    
	Elift=134.0;   
	furgc=17.0;     
	fusc=47.0;

% '>>>>>>>>Reactor Riser:'

	Mcpeff=10000.0; 
	Tbase=1100.0;   
	Tref=999.0;
	Tbasef=700.0;   
	dTstrp=35.0;    
	Hcoke=0.075;
	cpsv=0.80;      
	cpfv=0.81*(API/APInominal)^(-0.1);      
	cpfl=0.82*(API/APInominal)^(-0.1);   
	Qfr=309.0*(API/APInominal)^(-0.15);      
	Qsr=412.0;      
	rsarea=9.6;   
	hris=60.0;              
	rhov=0.57; 
    
% '>>>>>>>>Reactor Pressure:'

	k12=0.5; 
	dPfrac=9.5;	

%>>>>>>>>Define Variable for Common block
	xfcc=x;
	T3=xfcc(1); 			%Furnace firebox temperature
	T2=xfcc(2); 			%Temperature of fresh feed entering the reactor (Tpre IN MANUSCRIPT)

%(NEW) REACTOR WITH CONSTANT PRESSURE BECAUSE FRACTIONATOR SUBMODEL
%IS NOT CONSIDERED.
	
%	P7=xfcc(3); 			%WGC suction pressure
%	P5=xfcc(4); 			%Fractionator pressure (Pfra IN MANUSCRIPT)

% (NEW) THE FOLLOWING STATES WERE RECOUNTED
	P3=xfcc(3); 			%Lift air blower discharge pressure
	P6=xfcc(4); 			%Regenerator pressure (Preg IN MANUSCRIPT)
	rho=xfcc(5); 			%Density of catalyst in lift pipe
	P2=xfcc(6); 			%CAB discharge pressure
	Csc=xfcc(7); 			%Weight fraction of coke on spent catalyst
    Crgc=xfcc(8); 			%Weight fraction of coke on regenerated catalyst
    Treg=xfcc(9); 			%Temperature of regenerator (Treg IN MANUSCRIPT)
	Wsp=xfcc(10); 			%Inventory catalyst in regenerator standpipe
	Wreg=xfcc(11); 			%Inventory catalyst in regenerator
	Rn=xfcc(12); 			%Quantity of gas
	Wr=xfcc(13); 			%Inventory catalyst in reactor (Lrea IN MANUSCRIPT)
	Tr=xfcc(14); 			%Temperature of reactor riser (Trea IN MANUSCRIPT)
    Fair=xfcc(15); 			%Flow of air into regenerator 
	Pblp=xfcc(16); 			%Pressure at bottom of lift pipe
	P1=xfcc(17); 			%CAB suction pressure

% (NEW) CONTROL LOOPS WERE OPENED FOR THIS S-FUNCTION FOR WHICH
% THIS STATES NO LONGER EXIST
%	 PreHeatE=xfcc(20); 	%Preheater temperature controller integral error	
%    MainFracE=xfcc(21); 	%Fractionator pressure controller integral error
%    RegenE=xfcc(22); 		%Regenerator pressure controller integral error
%    RegenTE=xfcc(23); 		%Regenerator temperature controller integral error
%    CatIE=xfcc(24); 		%Reactor catalyst inventory controller integral error
%    ReacTE=xfcc(25); 		%Reactor temperature controller integral error

% (NEW) THE FOLLOWING STATES WERE RECOUNTED
    GFpvgo=xfcc(18);  		%CST MB for VGO
    GFp1=xfcc(19); 			%CST MB for PC1
    GFp2=xfcc(20); 			%CST MB for PC2
    GFp3=xfcc(21); 			%CST MB for PC3
    GFp4=xfcc(22); 			%CST MB for PC4
    GFc5=xfcc(23); 			%CST MB for C5+
    GFb=xfcc(24); 			%CST MB for Butane
    GFp=xfcc(25); 			%CST MB for Propane
    GFe=xfcc(26); 			%CST MB for Ethane
    GFm=xfcc(27); 			%CST MB for Methane

%(NEW) NO NEED FOR Fwg FILTER BECAUSE WET GAS COMPRESSOR IS NOT
%CONSIDERED.
%    Fwg=xfcc(28); 			%Filter for LPG

%(NEW) THE FOLLOWING STATES WERE RECOUNTED
    Fcoke=xfcc(28); 		%Filter for Coke

%(NEW) REDEFINING cont VARIABLE AS FUNCTION OF t
	cont=t/60;
	
   GTotal=[GFpvgo;GFp1;GFp2;GFp3;GFp4;GFc5;GFb;GFp;GFe;GFm];
   
%(NEW) REACTOR WITH CONSTANT PRESSURE BECAUSE FRACTIONATOR SUBMODEL
%IS NOT CONSIDERED.
   P4		= 24.9000000145994 + dPfrac;
   deltP	= P6 - P4;
   
 %>>>>>>>>>Define manipulated variables 
	ufcc=u(1:14);
	F3=ufcc(1); %Flow of fresh feed
	F4=ufcc(2); %Flow of slurry
	V12=ufcc(3); %Combustion air blower suction valve position
	V13=ufcc(4); %Combustion air blower vent valve position
	V14=ufcc(5); %Lift air blower vent valve position
	V15=ufcc(6); %Spill air valve position
	Vlift=ufcc(7); %Lift air blower steam
	V5=ufcc(8); %Wet gas flare valve position
	V16=ufcc(9); %WGC vent valve position 

% (NEW) CONTROL LOOPS WERE OPENED FOR THIS S-FUNCTION FOR WHICH
% THIS INPUTS NO LONGER EXIST
%   SPPH=ufcc(10); %Set point preheater (Temperature)
%   SPMF=ufcc(11); %Set point fractionator (Pressure) 
%   SPRG=ufcc(12); %Set point regenerator (Pressure) 
%   SPRGT=ufcc(13); %Set point regenerator (Temperature) 
%   SPCI=ufcc(14); %Set point reactor (Catalyst inventory)
%	SPRT=ufcc(15); %Set point reactor (Riser temperature) 

% (NEW) THE FOLLOWING INPUTS WERE RECOUNTED
% (NEW) THE DIRECT INPUT VARIABLES ARE NOW CONSIDERED AS INPUTS INSTEAD OF SET-POINTS
    p1a=ufcc(10); %Valve stem position valve for TC1 preheater (Temperature)

%(NEW) REACTOR WITH CONSTANT PRESSURE BECAUSE FRACTIONATOR SUBMODEL
%IS NOT CONSIDERED.
%    V4a=ufcc(11); %Valve stem position for PC2 fractionator (Pressure)

% (NEW) THE FOLLOWING INPUTS WERE RECOUNTED
    V7a=ufcc(11); %Valve stem position for PC1 regenerator (Pressure)
    p6a=ufcc(12); %Valve stem position for TC3 regenerator (Temperature)
    V3a=ufcc(13); %Valve stem position for LC1 reactor (Catalyst inventory)
	V2a=ufcc(14); %Valve stem position for TC2 reactor (Riser temperature)
    
    V13	=max(0.0d0,min(1.0d0,V13));
	V14	=max(0.0d0,min(1.0d0,V14));
	Vlift	=max(0.0d0,min(1.2d0,Vlift));

%>>>>>>>>>Nonlinear Valve Macro.'
 
	if(V12<=0.5)
   	V12=0.3*V12;
	else
   	V12=exp(-3.79424*(1.0d0-V12));
	end
  
	if(V13<=0.5)
 	 	V13=0.3*V13;
	else
  		V13=exp(-3.79424*(1.0d0-V13));
	end	

%*************************************************************'
%                  FCC UNIT                                   '
%*************************************************************'

%>>>>>>>>Independent variables

Tsc    = Tr - dTstrp;
Prgb   = P6 + Wreg/(144.0*rgarea);
 
%%%%%%%%%%%%%%%%% Controller Fractionator Pressure (PC2) %%%%%%%%%%%%%%%%%%
%(NEW) REACTOR WITH CONSTANT PRESSURE BECAUSE FRACTIONATOR SUBMODEL
%IS NOT CONSIDERED.
%alpha=1, hence, p(t)=valve(t)
%ePC2=(SPMF-P5)*-1;    

%KcPC2=1500;
%TaoPC2=200;                                                                         
%V4nom=50;

%if cont>=1   
%V4a=V4nom+KcPC2*(ePC2+MainFracE/TaoPC2);
%V4b=max(5,V4a);
%V4=min(95,V4b);
%else
%ePC2=0;    
%V4=V4nom; 
%end

%FV11   = (k11/2)*(V4/V4nom)^2*sqrt(max(0.0d0,P5 - P7)); 
%FV12   = k12*V5*sqrt(max(0.0d0,P5 - Patm));
%FV13   = k13*V16*Pvru;
%dMainFracE=ePC2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

%%%%%%%%%%%%%%%%% Controller Regenerator Temperature (TC3) %%%%%%%%%%%%%%%%
%(NEW) ORIGINAL PI LOOP WAS DEACTIVATED
%alpha=(100/160)
%eTC3=(SPRGT-Treg)*1;
    
%KcTC3=1.5;                    
%TaoTC3=800;            
p6nom=50;

%if cont>=1
%p6a=p6nom+KcTC3*(eTC3+RegenTE/TaoTC3);
p6b=max(8,p6a);
p6=min(152,p6b);
%else
%eTC3=0;    
%p6=p6nom;
%end

F7     = kca*(p6/p6nom)^1*sqrt(max(0.0d0,P2-Prgb)); %(Fa IN MANUSCRIPT)
V6=p6*(100/160);
%dRegenTE=eTC3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FV7    = k7*V13*sqrt(max(0.0d0,P2-Patm));
F10    = k9*V15*sqrt(max(0.0d0,P3-Prgb));

%(NEW) WET GAS COMPRESSOR NOT CONSIDERED
%Fwgi    =Flpg;

%%%%%%%%%%%%%%%%%%%%% Controller Regenerator Pressure (PC1)%%%%%%%%%%%%%%%%
%(NEW) ORIGINAL PI LOOP WAS DEACTIVATED
%alpha=1, hence, p(t)=valve(t)
%ePC1=(SPRG-P6)*-1;     

%KcPC1=1500;                             
%TaoPC1=200;                              
V7nom=50;

%if cont>=1
%V7a=V7nom+KcPC1*(ePC1+RegenE/TaoPC1);
V7b=max(10,V7a);
V7=min(95,V7b);
%else
%ePC1=0;    
%V7=V7nom;    
%end 

Fsg    = (1/2)*(V7/V7nom)*k14*sqrt(max(0.0d0,P6-Patm)); %(Ffg IN MANUSCRIPT)
%dRegenE=ePC1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhog   = 0.0933*P6/(Treg + Tf);   
F9     = klift*sqrt(max(0.0d0,P3 - Pblp));
vs     = ((Fsg + Fair)/2.0)*(1.0/(rhog*rgarea));
epsf   = 0.332 + 0.06*vs;
rhocdn = rhoprt*(1.0 - epsf);
rhocdl = -0.878 + 0.582*vs;
zbed   = min(zcyc,(2.85 + 0.8*vs + (Wreg - rhocdl*rgarea*...
               zcyc)/(rgarea*rhocdn))*(1.0/(1.0 - rhocdl/rhocdn)));
epse   = max(epsf,epsf + (1.904 + 0.363*vs - 0.048*vs*vs)/zbed);
epse   = min(1.0d0,epse);

%>>>>>>>Spent Catalyst Ubend.'

dPsc   = 144.0*(P4-Pblp) + ((Wr/Astrp) + (Estrp-Elift)*rhoc);

%%%%%%%%%%%%%%%%%% Controller Catalyst Inventory (LC1)%%%%%%%%%%%%%%%%%%%%%                   
%(NEW) ORIGINAL PI LOOP WAS DEACTIVATED
%alpha=1, hence, p(t)=valve(t)
%eLC1=(SPCI-Wr)*-1;

%KcLC1=0.0001;                       
%TaoLC1=1000;  
V3nom=50;

%if cont>=1
%V3a = V3nom+KcLC1*(eLC1+CatIE/TaoLC1);
V3b=max(5,V3a);
V3=min(95,V3b);
%else
%eLC1=0;    
%V3=V3nom;
%end

vsc    = (dPsc*areaus)/(Lus*fusc);
Fsc    = (V3/V3nom)^2*vsc*areaus*rhoc; %(Fsc IN MANUSCRIPT)
%dCatIE=eLC1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FH     = Fsc*(Csc - Crgc)*Hcoke; 

%>>>>>>>>Preheat System.'

DTin     = T3 - T1;
DTout    = T3 - T2;
Tlm      = (DTin - DTout)/(log(DTin/DTout));

%%%%%%%%%%%%%%%%%%%%%%% Controller Preheater (TC1)%%%%%%%%%%%%%%%%%%%%%%%%%
%(NEW) ORIGINAL PI LOOP WAS DEACTIVATED
%alpha=(118.75/100)
%eTC1=(SPPH-T2)*1;

%KcTC1=5;                          
%TaoTC1=100;                        
p1nom=50;

%if cont>=1   
%p1a=p1nom+KcTC1*(eTC1+PreHeatE/TaoTC1);
p1b=max(20,p1a); 
p1=min(80,p1b); 
%else
%eTC1=0;    
%p1=p1nom;
%end    
F5=(p1/p1nom)*F5nom; %(Ff IN MANUSCRIPT)
V1=p1*(118.75/100);
Qloss    = a1*F5*T3 - a2;
T2ss     = T1 + (UAf*Tlm/F3);
%dPreHeatE=eTC1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%>>>>>>>>>Wet Gas Compressor.'
%(NEW) WET GAS COMPRESSOR NOT CONSIDERED
%Crw      = Pvru/P7;
%Hwg      = 6884936.03341729*(Crw^2 - 1.0d0); 
%Fsucwg   = 106859.958742798 + sqrt(max(0.0d0,592504703729978-0.15*Hwg^2.0)); 
%F11      = 2.636d-6*Fsucwg*P7;

%>>>>>>>>>Lift Air Blower.'

M      = (kavg-1.0)/kavg/etapla;
sa     = samin + 1100.0*Vlift;
xp5p   = max(0.0d0,P3);     
Pbdla  = ((xp5p^M-Patm^M)*(sb/sa)^2.0 + Patm^M)^(1.0d0/M);
Fbsla  = 8600.0 + sqrt(max(0.0d0,2.582d8-1.068d5*Pbdla^2.0));
Fsucla = Fbsla*sa/sb; 
F8     = 0.0451*Patm*Fsucla/(Tatm+Tf);
FV8    = k8*V14*sqrt(max(0.0d0,P3-Patm));
rhoag  = 29.0*P6/R/(Tsc+Tf);
vairl  = F9/Alp/rhoag;
vcatl  = max(vairl-vslip,Fsc/Alp/rhoprt);

%>>>>>>>>>Combustion Air Blower.'

Pbca   = 14.7*P2/P1;
Fsucca = 45000.0+sqrt(max(0.0d0,1.5813d9-1.2491d6*Pbca^2.0));
F6     = 0.045*Fsucca*P1/(Tatm + Tf);
Fv6    = k6*V12*sqrt(Patm-P1);

%>>>>>>>Regenerated Catalyst Ubend.'
%%%%%%%%%%%%%%%%%%% Controller Reactor Temperature (TC2) %%%%%%%%%%%%%%%%%%
%(NEW) ORIGINAL PI LOOP WAS DEACTIVATED
%alpha=1, hence, p(t)=valve(t)
%eTC2=(SPRT-Tr)*1;

%KcTC2=0.5;                       
%TaoTC2=200;  
V2nom=50;

%if cont>=1
%V2a= V2nom+KcTC2*(eTC2+ReacTE/TaoTC2);
V2b=max(5,V2a);
V2=min(95,V2b);    
%else
%eTC2=0;    
%V2=V2nom;
%end  

eta1   = (Wsp/sparea) + (Etap-Eloi)*rhoc;
eta2   = areaur/(Lur*furgc);
eta3   = 144.0*eta2*(P6-P4) + eta2*eta1;
rho1   = rhov*rhoprt;
alpha1 = areaur*rhoc*eta3;
alpha2 = areaur*rhoc*eta2*hris;
beta1  = rhoprt*(F3+F4) + alpha2*rho1 - alpha1*rhov;
beta2  = (F3+F4)*(alpha2*rho1 - alpha1*rhoprt);
Frgc   = (V2/V2nom)^2*(-beta1 + ((beta1^2) - 4*rhov*beta2)^(0.5))/2.0/rhov; %(Frgc IN MANUSCRIPT)
%dReacTE=eTC2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vris   = (F3 + F4)/rhov + Frgc/rhoprt;
rhoris = (F3 + F4 + Frgc)/vris;
Wris   = (Frgc*rsarea*hris)/vris;
Prb    = P4 + (rhoris*hris)/144.0;

a0p=[2.24167820089635;0;0;0;0;0;0;0;0;0;0];%Initial Concentration of feed to reactor (mol/Kg gas)
[ypfr,taoF]=PFR(a0p,F3,F4,vris,API);%Riser model (PFR subroutine)

 dGFpvgo=taoF*(ypfr(1)-GFpvgo);% CST model
 dGFp1=taoF*(ypfr(2)-GFp1);% CST model
 dGFp2=taoF*(ypfr(3)-GFp2);% CST model
 dGFp3=taoF*(ypfr(4)-GFp3);% CST model
 dGFp4=taoF*(ypfr(5)-GFp4);% CST model
 dGFc5=taoF*(ypfr(6)-GFc5);% CST model
 dGFb=taoF*(ypfr(7)-GFb);% CST model
 dGFp=taoF*(ypfr(8)-GFp);% CST model
 dGFe=taoF*(ypfr(9)-GFe);% CST model
 dGFm=taoF*(ypfr(10)-GFm);% CST model

MWp5=446.094358949534;
MWp4=378.894679684485;
MWp3=292.185529267655;
MWp2=206.951432493898;
MWp1=120.941794467086;
MWc5=85.1347950005991;
MWb=58.12;
MWp=44.1;
MWe=30.07;
MWm=16.04;
MWc=400;

MWT=[MWp5;MWp4;MWp3;MWp2;MWp1;MWc5;MWb;MWp;MWe;MWm;MWc]; %(g/mol)
 
yfinal=([GTotal;ypfr(end)].*MWT)*(F3+F4)*(1/2.20462)*(1/453.592); %Mass of each component (lb/s)

Fcokei=yfinal(end);

%>>>>>>>>>Regenerator Calculations:'

zposition %Regenerator temperature profile and MB

Qair   = Fair*cpair*(Tair - Tbase);
Qh     = FH*delhH;
Qc     = Fair*(cocyc*delH1 + co2cyc*delH2);
Qsc    = Fsc*cpc*(Tsc-Tbase);
Qfg    = (Fair*(o2cyc*cpO2 + cocyc*cpCO + co2cyc*cpCO2 + 0.79*cpN2) + 0.5*FH*cpH2O)*(Tcyc - Tbase);
Qrgc   = Frgc*cpc*(Treg - Tbase);
Qin    = Qair + Qh + Qc + Qsc;
Qout   = Qfg + Qrgc + Qe;
hspp   = hsp-(Wsp/rhoc/sparea);
Fspadj = 20.0*(3.0 - min(3.0d0,hspp));
Fsp    = fof*(sparea^0.5)*(zbed - zsp) - 4925.0 - Fspadj;
splev  = Wsp/(rhoc*sparea);
Vregg  = rgarea*zcyc - rgarea*zbed*(1.0-epse);

dn     = Fair - Fsg;
dWc    = (Fsc*Csc-FH) - (Fsp*Crgc+12.0*Fair*(cocyc+co2cyc));
dTreg  = (Qin - Qout)/((Wreg + Wsp)*cpc + MI); 

%>>>>>>>Reactor Riser Energy Calculations'

dHcrak = 172.7  + 3.0*(Tr - Tref);
Qrcat  = Frgc*cpc*(Tr-Tbase);
Qslury = F4*(cpsv*(Tr-Tref)+Qsr);
Qff    = F3*(cpfv*(Tr-Tref)+Qfr);
Qcrak  = (F3 + F4)*dHcrak;
Qrin   = Qrgc + F3*cpfl*(T2-Tbasef);
Qrout  = Qrcat + Qslury + Qcrak + Qff;


dWr    = Frgc - Fsc;
dTr    = (Qrin - Qrout)/(Mcpeff);
dCsc   = (Frgc*Crgc + Fcoke - Fsc*Csc - Csc*dWr)/Wr;

%(NEW) REACTOR WITH CONSTANT PRESSURE BECAUSE FRACTIONATOR SUBMODEL
%IS NOT CONSIDERED.
%dP5    = 0.833*(Fwg - FV11 - FV12 + FV13);

%>>>>>>>Power and electric motor AMP Calculations'
FlowCAB=(((F6/18.01528)*(10.731557)*(Tatm+459.67))/P1)*60;%ft^3/min
PCAB=(((Tatm+459.67)/0.9)*(14.7/491.67)*FlowCAB*12^2*((P2/P1) - 1))/33000;%HP
ACAB=PCAB*746/(13800*0.9*1.1*1.73);%Amp

%(NEW) WET GAS COMPRESSOR NOT CONSIDERED
%FlowWGC=(((Fwg/453.59)*(10.731557)*(560.96))/P7)*60;%ft^3/min
%PWGC=((1/2)*(((Tcondenser-273.15)*1.8+491.67)/0.9)*(14.7/491.67)*FlowWGC*12^2*(Crw^2 - 1.0d0))/33000;%HP
%AWGC=PWGC*746/(13200*0.9*0.98*1.73);%Amp

%>>>>>>Define the output'
fac=xco*28+xco2*44+xo2*32+22.12;

%(NEW) REACTOR WITH CONSTANT PRESSURE BECAUSE FRACTIONATOR SUBMODEL
%IS NOT CONSIDERED.
%yp(1)=P4;

%(NEW) OUTPUTS RECOUNTED
yp(1)=deltP;
yp(2)=Fair*29.0d0;
yp(3)=P6;
yp(4)=T3;
yp(5)=T2;
yp(6)=Tr;
yp(7)=Treg;
yp(8)=splev;
yp(9)=Tcyc;
yp(10)=Tcyc-Treg;
yp(11)=xco*28*1.e+6/fac;
yp(12)=xo2*Fair*100/Fsg;
yp(13)=Csc;
yp(14)=Crgc;
yp(15)=Fair;
yp(16)=Wris;
yp(17)=Wreg;
yp(18)=Wsp;

%(NEW) FRACTIONATOR PRESSURE NOT CONSIDERED BECAUSE FRACTIONATOR 
%SUBMODEL AND WET GAS COMPRESSOR are NOT CONSIDERED.
%yp(20)=P5;
%yp(19)=V4;

%(NEW) OUTPUTS RECOUNTED
yp(19)=V6;
yp(20)=V7;
yp(21)=V3;
yp(22)=V1;
yp(23)=V2;
yp(24)=T2;

%(NEW) FRACTIONATOR PRESSURE NOT CONSIDERED BECAUSE FRACTIONATOR 
%SUBMODEL IS NOT CONSIDERED.
%yp(28)=P5;

%(NEW) OUTPUTS RECOUNTED
yp(25)=P6;
yp(26)=Tr;
yp(27)=Treg;
yp(28)=Wr;
yp(29)=Frgc*60;
yp(30)=Fsc*60;
yp(31)=Fcoke*60;
yp(32)=yfinal(1);
yp(33)=yfinal(2);
yp(34)=yfinal(3);
yp(35)=yfinal(4);
yp(36)=yfinal(5);
yp(37)=yfinal(6);
yp(38)=yfinal(7);
yp(39)=yfinal(8);
yp(40)=yfinal(9);
yp(41)=yfinal(10);
yp(42)=ones(1,11)*yfinal-F3;%Reactor MB
yp(43)=ACAB;%Power CAB

%(NEW) WET GAS COMPRESSOR NOT CONSIDERED
%yp(48)=AWGC;%Power WGC

%(NEW) OUTPUTS RECOUNTED
yp(44)=F5*60;
yp(45)=F7*60;
yp(46)=Fsg*60;

%(NEW) WET GAS COMPRESSOR NOT CONSIDERED
%yp(47)=FV11*60;

%*************************************************************'
%                  TIME DERIVATIVES                           '
%*************************************************************'

%>>>>>>>derivatives'
delta(1) = (F5*DHfu-UAf*Tlm-Qloss)/taufb;
delta(2) = (T2ss-T2)/taufo;

%(NEW) REACTOR WITH CONSTANT PRESSURE BECAUSE FRACTIONATOR SUBMODEL
%IS NOT CONSIDERED.
%delta(3) = 5.0*(FV11 - F11);
%delta(4) = dP5;

% (NEW) THE FOLLOWING STATES WERE RECOUNTED
delta(3) = (R*(Tdla+Tf)/29.0/Vdla)*(F8-FV8-F9-F10);
delta(4) = R*(Rn*dTreg+(Treg+Tf)*dn)/Vregg;
delta(5) = ((Fsc/vcatl/Alp)+rhoag-rho)/taufil;
delta(6) = R*(Tdca+Tf)*(F6-FV7-F7)/29.0/Vdca;
delta(7) = dCsc;
delta(8)= (dWc - Crgc*(Fsc-Fsp))/Wreg;
delta(9)= dTreg;
delta(10)= Fsp - Frgc;
delta(11)= Fsc - Fsp;
delta(12)= dn;
delta(13)= dWr;
delta(14)= dTr ;      
delta(15)= (-Fair+(F7+F9+F10)/29.0);
delta(16)= -Pblp+P6+(rho*hlift/144.0)+(zbed-zlp)*rhocdn/144.0;
delta(17)= R*(Tatm+Tf)*(Fv6-F6)/29.0/Vcms;

%(NEW) CONTROL LOOPS WERE OPENED FOR THIS S-FUNCTION FOR WHICH
% THIS STATES NO LONGER EXIST
%delta(20)=dPreHeatE;
%delta(21)=dMainFracE;
%delta(22)=dRegenE;
%delta(23)=dRegenTE;
%delta(24)=dCatIE;
%delta(25)=dReacTE;

% (NEW) THE FOLLOWING STATES WERE RECOUNTED
delta(18)=dGFpvgo;
delta(19)=dGFp1;
delta(20)=dGFp2;
delta(21)=dGFp3;
delta(22)=dGFp4;
delta(23)=dGFc5;
delta(24)=dGFb;
delta(25)=dGFp;
delta(26)=dGFe;
delta(27)=dGFm;

%(NEW) NO NEED FOR Fwg FILTER BECAUSE WET GAS COMPRESSOR IS NOT
%CONSIDERED.
%delta(28)=((Fwgi-Fwg))/5;

% (NEW) THE FOLLOWING STATES WERE RECOUNTED
delta(28)=(Fcokei-Fcoke)/5;

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS AND INITIAL CONDITIONS --------
    case 0,
sys(1)=28;%Number of continuous states
sys(2)=0;%Number of discrete states
sys(3)=46;%Number of outputs
sys(4)=14;%Number of inputs
sys(5)=0;%Reserved for root finding. Must be zero.
sys(6)=1;%Direct Feedthrough S-function parameter
sys(7)=1;%Number of sample times (at least one sample time is needed)

%----- Initial conditions:
xfcc=[ ...					%Old	New 	Description
	1564.06758445685;    ...%(1) 	(1)		Furnace firebox temperature
	616;                 ...%(2) 	(2) 	Temperature of fresh feed entering the reactor (Tpre IN MANUSCRIPT)
	...24.0440758965397; ...%(3) 	 		WGC suction pressure
	...24.9000000145994; ...%(4) 	 		Fractionator pressure (Pfra IN MANUSCRIPT)
	40.0420669348231;    ...%(5) 	(3) 	Lift air blower discharge pressure
	28.0000000041903;    ...%(6) 	(4) 	Regenerator pressure (Preg IN MANUSCRIPT)
	3.16075658859897;	 ...%(7) 	(5) 	Density of catalyst in lift pipe
	35.0875282505755;	 ...%(8) 	(6) 	CAB discharge pressure
	0.0101439627122109;	 ...%(9) 	(7) 	Weight fraction of coke on spent catalyst
	0.00250580481493606; ...%(10) 	(8)		Weight fraction of coke on regenerated catalyst
	1249.99997383664;    ...%(11) 	(9)		Temperature of regenerator (Treg IN MANUSCRIPT)
	9350.79072581430;	 ...%(12) 	(10)	Inventory catalyst in regenerator standpipe
	271210.435908199;	 ...%(13) 	(11)	Inventory catalyst in regenerator
	246.810877809621;	 ...%(14) 	(12)	Quantity of gas	
	98098.9991960623;	 ...%(15) 	(13)	Inventory catalyst in reactor (Lrea IN MANUSCRIPT)
	968.999859177516;    ...%(16) 	(14)	Temperature of reactor riser (Trea IN MANUSCRIPT)
	2.68045714195376;	 ...%(17) 	(15)	Flow of air into regenerator 
	30.4463683645626;    ...%(18) 	(16)	Pressure at bottom of lift pipe
	14.6380093682365;	 ...%(19) 	(17)	CAB suction pressure
	...-24.2765676899209;...%(20) 			Preheater temperature controller integral error	
	...-0.612986761842691;...%(21) 			Fractionator temperature controller integral error
	...0.611929610398922; ...%(22) 			Regenerator pressure controller integral error
	...-5641.47052429781; ...%(23) 			Regenerator temperature controller integral error
	...-30153150.1669079; ...%(24) 			Reactor catalyst inventory controller integral error
	...-1879.52204006479; ...%(25) 			Reactor temperature controller integral error
	0.00468987379029097; ...%(26)	(18) 	CST MB for VGO
	0.0481251434795964;  ...%(27) 	(19)	CST MB for PC1
	0.527273962836787;   ...%(28) 	(20)	CST MB for PC2
	0.357213788281919;   ...%(29) 	(21)	CST MB for PC3
	2.47339657486076;    ...%(30) 	(22)	CST MB for PC4
	1.41455228053760;    ...%(31) 	(23)	CST MB for C5+
	2.70179775410060;    ...%(32) 	(24)	CST MB for Butane
	2.44381704476869;    ...%(33) 	(25)	CST MB for Propane
	0.631637205219422;   ...%(34) 	(26)	CST MB for Ethane
	0.627951527686494;   ...%(35) 	(27)	CST MB for Methane
	...457.718461235902; ...%(36) 			Filter for LPG
	6.31203102548586];   ...%(37) 	(28)	Filter for Coke
	
x0=xfcc(:);

%----- Other building parameters:
str=[];     %str is always an empty matrix
ts=[0 0];   %Sample time array [0 0] for continuous system
simStateCompliance='UnknownSimState';
%--------------------------------------------------------------------
%----- STATE DERIVATIVES --------------------------------------------
    case 1,
sys=delta(:);
%--------------------------------------------------------------------
    case 2,
sys=[];
%----- OUPUTS -------------------------------------------------------
    case 3,
sys=yp(:);
%----- S-FUNCTION OTHER FLAGS ---------------------------------------
    case 4,
sys=[];
    case 9,
sys=[];
    otherwise   % Unexpected flags
DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end
%--------------------------------------------------------------------    
end
