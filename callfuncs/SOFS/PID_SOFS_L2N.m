function [KP,KI,KD,TI,TD,alp,P_,Sgm_FP,epsl_FP,Sgm_Pprp,epsl_Pprp]=...
	PID_SOFS_L2N(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,law)
%	Semidefinite programming (SDP) problem to obtain feedback PID 
%	controller gain which maximizes system decay rate with upper
%	bound on Hinf norm by L2 function norms. SDP seeks to 
%	maximize system decay rate solving the generalized eigenvalue 
%	problem (GEVP) using the gevp command.
%
%	The PID tuning is cast as a static output feedback (SOF) stabili-
%	zation problem. P_ is expressed as function of the control gain 
%	matrices, thus it is not a decision variable. After elimination of 
%	P_ as a decision, quadratic matrix terms in KP and KD arises which
%	are relaxed by an S-procedure.
%
% 	[KP,KI,KD,TI,TD,alp,P_,Sgm_FP,epsl_FP,Sgm_Pprp,epsl_Pprp]=...
%		PID_SOFS_L2N(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,law)
%
%---- Outputs -------------------------------------------------------
%   KP: 		Proportional action gain matrix
%   KI: 		Integral action gain matrix
%   KD: 		Derivative action gain matrix
%   TI: 		Integral time matrix
%   TD: 		Derivative time matrix
%   alp: 		Decay rate
%   P_: 		Augmented system Lyapunov matrix
%   Sgm_FP: 	Slack variable for S-procedure KP quad. term relax
%   epsl_FP:	Eigenvalue variable for S-procedure KP quad. term relax
%   Sgm_Pprp: 	Slack variable for S-procedure KD quad. term relax
%   epsl_Pprp: 	Eigenvalue variable for S-procedure KD quad. term relax
%--------------------------------------------------------------------
%---- Inputs --------------------------------------------------------
%	Poly: 	Structure storing LDI (A,Bu,Bw,Cy) polytope verteces
%		.Nx:		Number of states
%		.Ny:		Number of outputs
%		.Nw:		Number of exogenous inputs
%		.Nl:		Number of polytope vertices
%		.A{1:Nl}:	Autonomous system matrix (Nx,Nx)
%   	.Bu{1:Nl}: 	Maniputaled inputs gain matrix (Nx,Ny)
%   	.Bw{1:Nl}: 	Exogenous inputs gain matrix (Nx,Nw)
%   	.Cy{1:Nl}: 	Output-state matrix relationship (Ny,Nx)
%	gL2:	Upper bound on L2 norm
%   Tk: 	Control action structure matrix
%	PhiD:	Derivative kick preventor constant matrix (diagonal)
%   KPlb: 	Proportional gain matrix lower bound
%   KPub: 	Proportional gain matrix upper bound
%   TIlb: 	Integral time matrix lower bound
%   TIub: 	Integral time matrix upper bound
%   TDlb: 	Derivative time matrix lower bound
%   TDub: 	Derivative time matrix upper bound
%   law: 	PID law (default values: 'P', 'PI', 'PID')
%--------------------------------------------------------------------

%	Developed for the paper:
%	Trica, D. J. 2024. Multiple PID tuning strategy based on 
%	non-iterative linear matrix inequalities to solve large systems.
%
%   Diego Trica (diegotrica@gmail.com)
%--------------------------------------------------------------------
warning('off','MATLAB:rankDeficientMatrix');
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');

%Initialization parameters ------------------------------------------
rtol=1e-6;		%SDP relative tolerance
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces
Tk=sign(Tk); 	%Control action structure normalization
Bu=Poly.Bu;		%System maniputaled inputs polytope
Cy=Poly.Cy;		%Controlled output-state relationship polytope
Tl=cell(Nl,1);	%Left-multiplier matrix of equivalence relation
Tr=cell(Nl,1);	%Right-multiplier matrix of equivalence relation
R=cell(Nl,1);	%Weighting matrix of equivalence relation

for i=1:Nl
    Bu{i}=Bu{i}*Tk; %Bu columns permutation by Tk
	R{i}=Bu{i}'*Bu{i};
	Tl{i}=[Bu{i} null(Bu{i}')];
	Tr{i}=[Cy{i};null(Cy{i})']';
end

%Storing into Poly
Poly.Bu=Bu; %Bu columns permutation by Tk
Poly.Tl=Tl; Poly.Tr=Tr; Poly.R=R;

%K_ variable structure
switch law
    case 'P'
sK_=diag(1:Ny);

    case 'PI'
sK_=[diag(1:Ny),diag(Ny+1:2*Ny)];

    case 'PID'
sK_=[diag(1:Ny),diag(Ny+1:2*Ny),diag(2*Ny+1:3*Ny)];
end
%--------------------------------------------------------------------
%Set LMI constraints ------------------------------------------------
switch law
    case 'P'
LMISYS=LMI_P(Poly,gL2,sK_,KPlb,KPub);
    
	case 'PI'
LMISYS=LMI_PI(Poly,gL2,sK_,KPlb,KPub,TIlb,TIub);

    case 'PID'
LMISYS=LMI_PID(Poly,gL2,sK_,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub);
end 

%Solve SDP problem
LMIOPT=[rtol,2000,-1,50,0]; %LMI solver options
[alp,z]=gevp(LMISYS,Nl,LMIOPT,[],[],-Inf);
%--------------------------------------------------------------------
%Output solution ----------------------------------------------------
if isempty(alp)
	warning('Unfeasible problem. Check or modify constraints');
	KP=[];
	KI=[]; TI=[];
	KD=[]; TD=[];
	alp=[];
	P_=[];
	Sgm_FP=[]; epsl_FP=[];
	Sgm_Pprp=[]; epsl_Pprp=[];
else
%Augmented control gain matrix
K_=dec2mat(LMISYS,z,1);

if Ny<Nx
	%Augmented proportional gain matrix
	Kprp=dec2mat(LMISYS,z,6);
	
	%Slack variable for quadratic form disguise
	Sgm_FP=dec2mat(LMISYS,z,7);
	epsl_FP=dec2mat(LMISYS,z,8);
else
	%Augmented proportional gain matrix
	Kprp=dec2mat(LMISYS,z,3);

	%Slack variable for quadratic form disguise
	Sgm_FP=dec2mat(LMISYS,z,4);
	epsl_FP=dec2mat(LMISYS,z,5);
end

KI=[]; TI=[];
KD=[]; TD=[];
Sgm_Pprp=[]; epsl_Pprp=[];

%KP, KI, KD and P_
switch law
    case 'P'
KP=K_(:,1:Ny);
Pprp=cell(Nl,1);
for i=1:Nl
	Pprp{i}=0.5*(Tl{i}*Kprp*Tr{i}' + (Tl{i}*Kprp*Tr{i}')');
end
P_=Pprp;

%KP unpermutation by Tk
KP=Tk*KP;

    case 'PI'
KP=K_(:,1:Ny);
KI=K_(:,1+Ny:2*Ny);
TI=KP/KI;
Pprp=cell(Nl,1); P_=cell(Nl,1);
for i=1:Nl
	Pprp{i}=0.5*(Tl{i}*Kprp*Tr{i}' + (Tl{i}*Kprp*Tr{i}')');
	P_{i}=[Pprp{i}    Bu{i}
		   Bu{i}'     eye(Ny)];
end

%KP, KI, TI unpermutation by Tk
KP=Tk*KP;
KI=Tk*KI; TI=abs(Tk)*TI;

    case 'PID'
if Ny<Nx
	%Slack variable for quadratic form disguise
	Sgm_Pprp=dec2mat(LMISYS,z,9);
	epsl_Pprp=dec2mat(LMISYS,z,10);
else
	%Slack variable for quadratic form disguise
	Sgm_Pprp=dec2mat(LMISYS,z,6);
	epsl_Pprp=dec2mat(LMISYS,z,7);	
end
		
KP=K_(:,1:Ny);
KI=K_(:,1+Ny:2*Ny);
TI=KP/KI;
KD=K_(:,1+2*Ny:3*Ny);
TD=KD/KP;
Pprp=cell(Nl,1); SD=cell(Nl,1); P_=cell(Nl,1);
for i=1:Nl
	SD{i}=Bu{i}*KD - Cy{i}'*(KD'*KD);

	Pprp{i}=(Tl{i}*Kprp*Tr{i}' + (Tl{i}*Kprp*Tr{i}')' ...
		- (Cy{i}'*PhiD*SD{i}') - (Cy{i}'*PhiD*SD{i}')')/2;

	P_{i}=[Pprp{i}           Bu{i}   	(Bu{i}-Cy{i}'*KD')
		   Bu{i}'         	 eye(Ny)     zeros(Ny)
		   (Bu{i}'-KD*Cy{i}) zeros(Ny)	 inv(PhiD)];
end

%KP, KI, TI, KD and TD unpermutation by Tk
KP=Tk*KP;
KI=Tk*KI; TI=abs(Tk)*TI;
KD=Tk*KD; TD=abs(Tk)*TD;

end
end
%--------------------------------------------------------------------
warning('on','MATLAB:rankDeficientMatrix');
warning('on','MATLAB:singularMatrix');
warning('on','MATLAB:nearlySingularMatrix');
end

function LMISYS=LMI_P(Poly,gL2,sK_,KPlb,KPub)
%Parameters
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces
A=Poly.A;		%Autonomous system matrix polytope
%Bu=Poly.Bu;	%System maniputaled inputs polytope
Bw=Poly.Bw;		%System exogenous inputs polytope
Cy=Poly.Cy;		%Controlled output-state relationship polytope
Tl=Poly.Tl;     %Left-multiplier matrix of equivalence relation
Tr=Poly.Tr;     %Right-multiplier matrix of equivalence relation
R=Poly.R;       %Weighting matrix of equivalence relation

%Initialize LMI system setting
setlmis([]);
	
%LMI system variables
lmivar(3,sK_);  					%Augmented control gain
[KP,~,sKP]=lmivar(3,sK_(:,1:Ny));	%Proportional gain

%Structure for augmented propotional gain
if Ny<Nx
	[~,~,sXP]=lmivar(2,[Ny,Nx-Ny]);
	[~,~,sYP]=lmivar(2,[Nx-Ny,Ny]);
	[~,~,sZP]=lmivar(2,[Nx-Ny,Nx-Ny]);
else
	sXP=[]; sYP=[]; sZP=[];
end
	
%Augmented prop. gain	
Kprp=lmivar(3,[[sKP,sXP];[sYP,sZP]]);

%Slack variable for quadratic form disguise
sSgm_FP=[];
for i=1:Ny
	sSgm_FP=[sSgm_FP;[1 0]];
end
Sgm_FP=lmivar(1,sSgm_FP);
epsl_FP=lmivar(1,[1 1]);

%LMI system constraints
neq=0;
for i=1:Nl
%P_{i} > 0
	neq=neq+1;
	lmiterm([-neq 1 1 Kprp],0.5*Tl{i},Tr{i}','s');

%SgmF_{i} > 0
	MP=KPub'*R{i}*KPub; sqrtMP=sqrtm(MP);
	
	neq=neq+1;
	lmiterm([-neq 1 1 epsl_FP],1,1);
	lmiterm([-neq 1 1 Sgm_FP],-1,1);
	lmiterm([-neq 2 1 KP],sqrtm(R{i}),eye(Ny)/sqrtMP);
	lmiterm([-neq 2 1 epsl_FP],-sqrtm(R{i}),KPub/sqrtMP);
	lmiterm([-neq 2 2 0],1);

%H{i} < 0
	neq=neq+1;
	lmiterm([neq 1 1 Kprp],0.5*A{i}'*Tl{i},Tr{i}','s');
	lmiterm([neq 1 1 Kprp],Tl{i},0.5*Tr{i}'*A{i},'s');
	lmiterm([neq 1 1 0],Cy{i}'*Cy{i});
	lmiterm([neq 1 1 Sgm_FP],-(sqrtMP*Cy{i})',sqrtMP*Cy{i},'s');
	lmiterm([neq 2 1 Kprp],0.5*Bw{i}'*Tl{i},Tr{i}');
	lmiterm([neq 2 1 -Kprp],0.5*Bw{i}'*Tr{i},Tl{i}');
	lmiterm([neq 2 2 0],-gL2^2);	
end
	
%|KP| > |KPlb|
neq=neq+1;
lmiterm([-neq 1 1 KP],1/2,1,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([neq 1 1 0],(KPlb+KPlb')/2); %Forcing symmetric form to not disp warning in command window
	
%|KP| < |KPub|
neq=neq+1;
lmiterm([neq 1 1 KP],1/2,1,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([-neq 1 1 0],(KPub+KPub')/2); %Forcing symmetric form to not disp warning in command window

%Sgm_FP > 0
neq=neq+1;
lmiterm([-neq 1 1 Sgm_FP],1,1);

%epsl_FP < 1
neq=neq+1;
lmiterm([neq 1 1 epsl_FP],1,1);
lmiterm([-neq 1 1 0],1);

for i=1:Nl
%(LFC) F{i} < 2*alp*P_{i}
	MP=KPub'*R{i}*KPub; sqrtMP=sqrtm(MP);

	neq=neq+1;
	lmiterm([neq 1 1 Kprp],0.5*A{i}'*Tl{i},Tr{i}','s');
	lmiterm([neq 1 1 Kprp],Tl{i},0.5*Tr{i}'*A{i},'s');
	lmiterm([neq 1 1 Sgm_FP],-(sqrtMP*Cy{i})',sqrtMP*Cy{i},'s');
	lmiterm([-neq 1 1 Kprp],Tl{i},Tr{i}','s');
end 
	
%End and store LMI system
LMISYS=getlmis;
end

function LMISYS=LMI_PI(Poly,gL2,sK_,KPlb,KPub,TIlb,TIub)
%Parameters
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces
A=Poly.A;		%Autonomous system matrix polytope
Bu=Poly.Bu;     %System maniputaled inputs polytope
Bw=Poly.Bw;		%System exogenous inputs polytope
Cy=Poly.Cy;     %Controlled output-state relationship polytope
Tl=Poly.Tl;     %Left-multiplier matrix of equivalence relation
Tr=Poly.Tr;     %Right-multiplier matrix of equivalence relation
R=Poly.R;       %Weighting matrix of equivalence relation

%Initialize LMI system setting
setlmis([]);
	
%LMI system variables
lmivar(3,sK_);  					%Augmented control gain
[KP,~,sKP]=lmivar(3,sK_(:,1:Ny));	%Proportional gain

%Structure for augmented propotional gain
if Ny<Nx
	[~,~,sXP]=lmivar(2,[Ny,Nx-Ny]);
	[~,~,sYP]=lmivar(2,[Nx-Ny,Ny]);
	[~,~,sZP]=lmivar(2,[Nx-Ny,Nx-Ny]);
else
	sXP=[]; sYP=[]; sZP=[];
end
	
%Augmented prop. gain	
Kprp=lmivar(3,[[sKP,sXP];[sYP,sZP]]);

%Slack variable for quadratic form disguise
sSgm_FP=[];
for i=1:Ny
	sSgm_FP=[sSgm_FP;[1 0]];
end
Sgm_FP=lmivar(1,sSgm_FP);
epsl_FP=lmivar(1,[1 1]);

KI=lmivar(3,sK_(:,1+Ny:2*Ny));		%Integral gain

%LMI system constraints
neq=0;
for i=1:Nl
%P_{i} > 0
	neq=neq+1;
	lmiterm([-neq 1 1 Kprp],0.5*Tl{i},Tr{i}','s');
	lmiterm([-neq 2 1 0],Bu{i}');
	lmiterm([-neq 2 2 0],eye(Ny));

%SgmF_{i} > 0
	MP=KPub'*R{i}*KPub; sqrtMP=sqrtm(MP);
	
	neq=neq+1;
	lmiterm([-neq 1 1 epsl_FP],1,1);
	lmiterm([-neq 1 1 Sgm_FP],-1,1);
	lmiterm([-neq 2 1 KP],sqrtm(R{i}),eye(Ny)/sqrtMP);
	lmiterm([-neq 2 1 epsl_FP],-sqrtm(R{i}),KPub/sqrtMP);
	lmiterm([-neq 2 2 0],1);

%H{i} < 0
	neq=neq+1;
	lmiterm([neq 1 1 Kprp],0.5*A{i}'*Tl{i},Tr{i}','s');
	lmiterm([neq 1 1 Kprp],Tl{i},0.5*Tr{i}'*A{i},'s');
	lmiterm([neq 1 1 KI],Bu{i},Cy{i},'s');		
	lmiterm([neq 1 1 0],Cy{i}'*Cy{i});
	lmiterm([neq 1 1 Sgm_FP],-(sqrtMP*Cy{i})',sqrtMP*Cy{i},'s');
	lmiterm([neq 2 1 0],Bu{i}'*A{i});
	lmiterm([neq 2 1 KP],-2*R{i},Cy{i});
	lmiterm([neq 2 1 KI],1,Cy{i});
	lmiterm([neq 2 2 0],-2*R{i});
	lmiterm([neq 3 1 Kprp],0.5*Bw{i}'*Tl{i},Tr{i}');
	lmiterm([neq 3 1 -Kprp],0.5*Bw{i}'*Tr{i},Tl{i}');
	lmiterm([neq 3 2 0],Bw{i}'*Bu{i});
	lmiterm([neq 3 3 0],-gL2^2);
end

%|KP| > |KPlb|
neq=neq+1;
lmiterm([-neq 1 1 KP],1/2,1,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([neq 1 1 0],(KPlb+KPlb')/2); %Forcing symmetric form to not disp warning in command window
	
%|KP| < |KPub|
neq=neq+1;
lmiterm([neq 1 1 KP],1/2,1,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([-neq 1 1 0],(KPub+KPub')/2); %Forcing symmetric form to not disp warning in command window

%|KI| > |TIub\KP|
neq=neq+1;
lmiterm([-neq 1 1 KI],1/2,1,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([neq 1 1 KP],eye(Ny)/TIub,1/2,'s'); %Forcing symmetric form to not disp warning in command window
	
%|KI| < |TIlb\KP|
neq=neq+1;
lmiterm([neq 1 1 KI],1/2,1,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([-neq 1 1 KP],eye(Ny)/TIlb,1/2,'s'); %Forcing symmetric form to not disp warning in command window

%Sgm_FP > 0
neq=neq+1;
lmiterm([-neq 1 1 Sgm_FP],1,1);

%epsl_FP < 1
neq=neq+1;
lmiterm([neq 1 1 epsl_FP],1,1);
lmiterm([-neq 1 1 0],1);

for i=1:Nl
%(LFC) F{i} < 2*alp*P_{i}
	neq=neq+1;
	lmiterm([neq 1 1 Kprp],0.5*A{i}'*Tl{i},Tr{i}','s');
	lmiterm([neq 1 1 Kprp],Tl{i},0.5*Tr{i}'*A{i},'s');
	lmiterm([neq 1 1 KI],Bu{i},Cy{i},'s');
	lmiterm([neq 1 1 Sgm_FP],-(sqrtMP*Cy{i})',sqrtMP*Cy{i},'s');
	lmiterm([neq 2 1 0],Bu{i}'*A{i});
	lmiterm([neq 2 1 KP],-2*R{i},Cy{i});
	lmiterm([neq 2 1 KI],1,Cy{i});
	lmiterm([neq 2 2 0],-2*R{i});
	lmiterm([-neq 1 1 Kprp],Tl{i},Tr{i}','s');
	lmiterm([-neq 2 1 0],2*Bu{i}');
	lmiterm([-neq 2 2 0],2*eye(Ny));
end

%End and store LMI system
LMISYS=getlmis;
end

function LMISYS=...
	LMI_PID(Poly,gL2,sK_,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub)
%Parameters
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces
A=Poly.A;		%Autonomous system matrix polytope
Bu=Poly.Bu;     %System maniputaled inputs polytope
Bw=Poly.Bw;		%System exogenous inputs polytope
Cy=Poly.Cy;     %Controlled output-state relationship polytope
Tl=Poly.Tl;     %Left-multiplier matrix of equivalence relation
Tr=Poly.Tr;     %Right-multiplier matrix of equivalence relation
R=Poly.R;       %Weighting matrix of equivalence relation

%Initialize LMI system setting
setlmis([]);
	
%LMI system variables
lmivar(3,sK_);  					%Augmented control gain
[KP,~,sKP]=lmivar(3,sK_(:,1:Ny));	%Proportional gain

if Ny<Nx
	[~,~,sXP]=lmivar(2,[Ny,Nx-Ny]);
	[~,~,sYP]=lmivar(2,[Nx-Ny,Ny]);
	[~,~,sZP]=lmivar(2,[Nx-Ny,Nx-Ny]);
else
	sXP=[]; sYP=[]; sZP=[];
end
	
%Augmented prop. gain	
Kprp=lmivar(3,[[sKP,sXP];[sYP,sZP]]);

%Slack variable for quadratic form disguise
sSgm_FP=[]; sSgm_Pprp=[];
for i=1:Ny
	sSgm_FP=[sSgm_FP;[1 0]];
	sSgm_Pprp=[sSgm_Pprp;[1 0]];
end
Sgm_FP=lmivar(1,sSgm_FP);
epsl_FP=lmivar(1,[1 1]);

MD=(TDub*KPub)'*PhiD*(TDub*KPub); sqrtMD=sqrtm(MD);
Sgm_Pprp=lmivar(1,sSgm_Pprp);
epsl_Pprp=lmivar(1,[1 1]);

KI=lmivar(3,sK_(:,1+Ny:2*Ny));		%Integral gain
KD=lmivar(3,sK_(:,1+2*Ny:3*Ny));	%Derivative gain

%LMI system constraints
neq=0;
for i=1:Nl
%P_{i} > 0
	neq=neq+1;
	lmiterm([-neq 1 1 Kprp],0.5*Tl{i},Tr{i}','s');
	lmiterm([-neq 1 1 KD],-0.5*Bu{i},PhiD*Cy{i},'s');
	lmiterm([-neq 1 1 Sgm_Pprp],0.5*(sqrtMD*Cy{i})',sqrtMD*Cy{i},'s');
	lmiterm([-neq 2 1 0],Bu{i}');
	lmiterm([-neq 3 1 0],Bu{i}');
	lmiterm([-neq 3 1 KD],-1,Cy{i});
	lmiterm([-neq 2 2 0],eye(Ny));
	lmiterm([-neq 3 3 0],eye(Ny)/PhiD);

%SgmF_{i} > 0
	MP=KPub'*R{i}*KPub; sqrtMP=sqrtm(MP);
	
	neq=neq+1;
	lmiterm([-neq 1 1 epsl_FP],1,1);
	lmiterm([-neq 1 1 Sgm_FP],-1,1);
	lmiterm([-neq 2 1 KP],sqrtm(R{i}),eye(Ny)/sqrtMP);
	lmiterm([-neq 2 1 epsl_FP],-sqrtm(R{i}),KPub/sqrtMP);
	lmiterm([-neq 2 2 0],1);

%H{i} < 0
	neq=neq+1;
	lmiterm([neq 1 1 Kprp],0.5*A{i}'*Tl{i},Tr{i}','s');
	lmiterm([neq 1 1 Kprp],Tl{i},0.5*Tr{i}'*A{i},'s');
	lmiterm([neq 1 1 KD],0.5*Bu{i},PhiD*Cy{i}*A{i},'s');
	lmiterm([neq 1 1 KD],-0.5*A{i}'*Bu{i},PhiD*Cy{i},'s');
	lmiterm([neq 1 1 KI],Bu{i},Cy{i},'s');		
	lmiterm([neq 1 1 0],Cy{i}'*Cy{i});
	lmiterm([neq 1 1 Sgm_FP],-(sqrtMP*Cy{i})',sqrtMP*Cy{i},'s');
	lmiterm([neq 2 1 0],Bu{i}'*A{i});
	lmiterm([neq 2 1 KP],-2*R{i},Cy{i});
	lmiterm([neq 2 1 KI],1,Cy{i});
	lmiterm([neq 3 1 0],Bu{i}'*A{i});
	lmiterm([neq 3 1 0],-PhiD*Bu{i}');
	lmiterm([neq 3 1 KP],-2*R{i},Cy{i});
	lmiterm([neq 3 1 KD],PhiD,Cy{i});
	lmiterm([neq 2 2 0],-2*R{i});
	lmiterm([neq 3 2 0],-2*R{i});
	lmiterm([neq 3 3 0],-2*(R{i}+eye(Ny)));
	lmiterm([neq 4 1 Kprp],0.5*Bw{i}'*Tl{i},Tr{i}');
	lmiterm([neq 4 1 -Kprp],0.5*Bw{i}'*Tr{i},Tl{i}');
	lmiterm([neq 4 1 KD],-0.5*Bw{i}'*Bu{i},PhiD*Cy{i});
	lmiterm([neq 4 1 -KD],0.5*Bw{i}'*Cy{i}'*PhiD,Bu{i}');
	lmiterm([neq 4 2 0],Bw{i}'*Bu{i});
	lmiterm([neq 4 3 0],Bw{i}'*Bu{i});
	lmiterm([neq 4 4 0],-gL2^2);
end

%|KP| > |KPlb|
neq=neq+1;
lmiterm([-neq 1 1 KP],1/2,1,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([neq 1 1 0],(KPlb+KPlb')/2); %Forcing symmetric form to not disp warning in command window
	
%|KP| < |KPub|
neq=neq+1;
lmiterm([neq 1 1 KP],1/2,1,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([-neq 1 1 0],(KPub+KPub')/2); %Forcing symmetric form to not disp warning in command window

%|KI| > |TIub\KP|
neq=neq+1;
lmiterm([-neq 1 1 KI],1/2,1,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([neq 1 1 KP],eye(Ny)/TIub,1/2,'s'); %Forcing symmetric form to not disp warning in command window
	
%|KI| < |TIlb\KP|
neq=neq+1;
lmiterm([neq 1 1 KI],1/2,1,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([-neq 1 1 KP],eye(Ny)/TIlb,1/2,'s'); %Forcing symmetric form to not disp warning in command window

%|KD| > |TDlb*KP|
neq=neq+1;
lmiterm([-neq 1 1 KD],1/2,1,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([neq 1 1 KP],TDlb,1/2,'s'); %Forcing symmetric form to not disp warning in command window
	
%|KD| < |TDub*KP|
neq=neq+1;
lmiterm([neq 1 1 KD],1/2,1,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([-neq 1 1 KP],TDub,1/2,'s'); %Forcing symmetric form to not disp warning in command window

%Sgm_Pprp > 0
neq=neq+1;
lmiterm([-neq 1 1 Sgm_Pprp],1,1);

%epsl_Pprp < 1
neq=neq+1;
lmiterm([neq 1 1 epsl_Pprp],1,1);
lmiterm([-neq 1 1 0],1);

%SgmPprp_ > 0
	neq=neq+1;
	lmiterm([-neq 1 1 epsl_Pprp],1,1);
	lmiterm([-neq 1 1 Sgm_Pprp],-1,1);
	lmiterm([-neq 2 1 KD],sqrtm(PhiD),eye(Ny)/sqrtMD);
	lmiterm([-neq 2 1 epsl_Pprp],-sqrtm(PhiD),(KPub*TDub)/sqrtMD);
	lmiterm([-neq 2 2 0],1);

%Sgm_FP > 0
neq=neq+1;
lmiterm([-neq 1 1 Sgm_FP],1,1);

%epsl_FP < 1
neq=neq+1;
lmiterm([neq 1 1 epsl_FP],1,1);
lmiterm([-neq 1 1 0],1);

for i=1:Nl
%(LFC) F{i} < 2*alp*P_{i}
	MP=KPub'*R{i}*KPub; sqrtMP=sqrtm(MP);

	neq=neq+1;
	lmiterm([neq 1 1 Kprp],0.5*A{i}'*Tl{i},Tr{i}','s');
	lmiterm([neq 1 1 Kprp],Tl{i},0.5*Tr{i}'*A{i},'s');
	lmiterm([neq 1 1 KD],0.5*Bu{i},PhiD*Cy{i}*A{i},'s');
	lmiterm([neq 1 1 KD],-0.5*A{i}'*Bu{i},PhiD*Cy{i},'s');
	lmiterm([neq 1 1 KI],Bu{i},Cy{i},'s');
	lmiterm([neq 1 1 Sgm_FP],-(sqrtMP*Cy{i})',sqrtMP*Cy{i},'s');
	lmiterm([neq 2 1 0],Bu{i}'*A{i});
	lmiterm([neq 2 1 KP],-2*R{i},Cy{i});
	lmiterm([neq 2 1 KI],1,Cy{i});
	lmiterm([neq 3 1 0],Bu{i}'*A{i});
	lmiterm([neq 3 1 0],-PhiD*Bu{i}');
	lmiterm([neq 3 1 KP],-2*R{i},Cy{i});
	lmiterm([neq 3 1 KD],PhiD,Cy{i});
	lmiterm([neq 2 2 0],-2*R{i});
	lmiterm([neq 3 2 0],-2*R{i});
	lmiterm([neq 3 3 0],-2*(R{i}+eye(Ny)));
	lmiterm([-neq 1 1 Kprp],Tl{i},Tr{i}','s');
	lmiterm([-neq 1 1 KD],-Bu{i},PhiD*Cy{i},'s');
	lmiterm([-neq 1 1 Sgm_Pprp],(sqrtMD*Cy{i})',sqrtMD*Cy{i},'s');
	lmiterm([-neq 2 1 0],2*Bu{i}');
	lmiterm([-neq 3 1 0],2*Bu{i}');
	lmiterm([-neq 3 1 KD],-2,Cy{i});
	lmiterm([-neq 2 2 0],2*eye(Ny));
	lmiterm([-neq 3 3 0],2*eye(Ny)/PhiD);	
end

%End and store LMI system
LMISYS=getlmis;
end