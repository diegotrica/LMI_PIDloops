function [KP,KI,KD,TI,TD,alp,P_,Sgm_FP,epsl_FP,Sgm_Pprp,epsl_Pprp]=...
	PID_SOFS_L2N_YALMIP(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,law)
%	Semidefinite programming (SDP) problem to obtain feedback PID 
%	controller gain which maximizes system decay rate with upper
%	bound on Hinf norm by L2 function norms using YALMIP interface.
%	SDP seeks to maximize system decay rate solving the generalized 
%	eigenvalue problem (GEVP) using YALMIP bisection built-in function.
%
%	The PID tuning is cast as a static output feedback (SOF) stabili-
%	zation problem. P_ is expressed as function of the control gain 
%	matrices, thus it is not a decision variable. After elimination of 
%	P_ as a decision, quadratic matrix terms in KP and KD arises which
%	are relaxed by an S-procedure. 
%
% 	[KP,KI,KD,TI,TD,alp,P_,Sgm_FP,epsl_FP,Sgm_Pprp,epsl_Pprp]=...
%		PID_SOFS_L2N_YALMIP(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,law)
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
yalmip('clear');
warning('off','MATLAB:rankDeficientMatrix');
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');

%Initialization parameters ------------------------------------------
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
Poly.Bu=Bu;
Poly.Tl=Tl; Poly.Tr=Tr; Poly.R=R;

%Solve LMI problem --------------------------------------------------
switch law
    case 'P'
[vars,lmis,obj]=LMI_P(Poly,gL2,KPlb,KPub);
    
	case 'PI'
[vars,lmis,obj]=LMI_PI(Poly,gL2,KPlb,KPub,TIlb,TIub);

    case 'PID'
[vars,lmis,obj]=LMI_PID(Poly,gL2,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub);
end 

%YALMIP SDP setting options
% opt=sdpsettings(...
%     'solver','bisection',...
%     'bisection.absgaptol',1e-6, ...
%     'bisection.relgaptol',1e-3,...
%     'bisection.solver','sedumi',...
% 	  'sedumi.eps',1e-8,...
%     'sedumi.cg.qprec',1,'sedumi.cg.maxiter',49,...
%     'sedumi.stepdif',2,...'sedumi.stepdif',1,...'sedumi.stepdif',0,...
%     'sedumi.free',1,...'sedumi.free',0,...
%     'sedumi.sdp',1,...'sedumi.sdp',0,...
%     'verbose',1);
opt=sdpsettings(...
    'solver','bisection',...
    'bisection.absgaptol',1e-5, ...
    'bisection.relgaptol',1e-3,...
    'bisection.solver','sdpt3',...
    'verbose',1);

%Repeat feasibility problem with fixed alp, obtaining optimal alp by
%bisection method
stat=optimize(lmis,obj,opt);

if stat.problem==1
	warning(stat.info);
end
%--------------------------------------------------------------------
%Output solution ----------------------------------------------------
if stat.problem==1
	KP=[];
	KI=[]; TI=[];
	KD=[]; TD=[];
	alp=[];
	P_=[];
	Sgm_FP=[]; epsl_FP=[];
    Sgm_Pprp=[]; epsl_Pprp=[];
else
    
	alp=value(vars.alp);
	Sgm_FP=value(vars.Sgm_FP);
	epsl_FP=value(vars.epsl_FP);
	
    KI=[]; TI=[];
    KD=[]; TD=[];
    Sgm_Pprp=[]; epsl_Pprp=[];
    
	switch law
        case 'P'
    KP=value(vars.KP);
    P_=cell(Nl,1);
    for i=1:Nl
        P_{i}=value(vars.P_{i});
    end
    
    %KP unpermutation by Tk
    KP=Tk*KP;
	
        case 'PI'
    KP=value(vars.KP);
    KI=value(vars.KI); TI=KP/KI;
    P_=cell(Nl,1);
    for i=1:Nl
        P_{i}=value(vars.P_{i});
    end
    
    %KP, KI, TI unpermutation by Tk
    KP=Tk*KP;
    KI=Tk*KI; TI=abs(Tk)*TI;
	
        case 'PID'
    KP=value(vars.KP);
    KI=value(vars.KI); TI=KP/KI;
    KD=value(vars.KD); TD=KD/KP;	
    epsl_Pprp=value(vars.epsl_Pprp);
    Sgm_Pprp=value(vars.Sgm_Pprp);

    Kprp=value(vars.Kprp);
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

function [vars,lmis,obj]=LMI_P(Poly,gL2,KPlb,KPub)
%Parameters ---------------------------------------------------------
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nw=Poly.Nw;		%Exogenous input variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces
A=Poly.A;		%Autonomous system matrix polytope
%Bu=Poly.Bu;	%System maniputaled inputs polytope
Bw=Poly.Bw;		%System exogenous inputs polytope
Cy=Poly.Cy;	    %Controlled output-state relationship polytope
Tl=Poly.Tl;     %Left-multiplier matrix of equivalence relation
Tr=Poly.Tr;     %Right-multiplier matrix of equivalence relation
R=Poly.R;      %Weighting matrix of equivalence relation
	
%SDP variables ------------------------------------------------------
KP=sdpvar(Ny,Ny,'diagonal');

%Augmented propotional gain
if Ny<Nx
	XP=sdpvar(Ny,Nx-Ny,'full');
	YP=sdpvar(Nx-Ny,Ny,'full');
	ZP=sdpvar(Nx-Ny,Nx-Ny,'full');
	Kprp=[[KP,XP];[YP,ZP]];
else
	Kprp=KP;
end

%Slack variable for quadratic form disguise
Sgm_FP=sdpvar(Ny,Ny,'diagonal');
epsl_FP=sdpvar(1);

for i=1:Nl
%Lyapunov matrix proportional term
	Gprp{i}=(Tl{i}*Kprp*Tr{i}') + (Tl{i}*Kprp*Tr{i}')';
	Pprp{i}=Gprp{i}/2;

%Entire Lyapunov matrix
	P_{i}=Pprp{i};
end

%Eigenvalue
alp=sdpvar(1);

%Storing SDP variables
vars.KP=KP;
vars.Kprp=Kprp; vars.Pprp=Pprp;
vars.P_=P_;
vars.alp=alp;
vars.Sgm_FP=Sgm_FP;
vars.epsl_FP=epsl_FP;

%LMI system constraints ---------------------------------------------
lmis=[...
	KP >= KPlb, KP <= KPub,...
	Sgm_FP >= eps*eye(Ny),...
	epsl_FP <= 1 ...	
	];

for i=1:Nl
% Slack variable LMI constraint
	MP=KPub'*R{i}*KPub; sqrtMP=sqrtm(MP);
    pencilP=(sqrtm(R{i})*(KP-epsl_FP*KPub)/sqrtMP);
    SgmF_{i}=[[epsl_FP*eye(Ny)-Sgm_FP, pencilP'];
              [pencilP               , eye(Ny) ]];

%F_(KP)
	FP=...	%FP with quadratic form disguise using slack variable Sgm_FP
        ( (A{i}'*Gprp{i}) + (A{i}'*Gprp{i})' )/2 ...
        -2*(sqrtMP*Cy{i})'*Sgm_FP*(sqrtMP*Cy{i});    

%H_(KP):    
    HP=FP + (Cy{i}'*Cy{i} + eps*ones(size(Cy{i}'*Cy{i},1))); %adding an almost null matrix to handle numerical issues
    SP=Bw{i}'*Gprp{i}/2;

    H_{i}=[[HP, SP'           ];...
           [SP, -gL2^2*eye(Nw)]];

    lmis=[lmis,...
	    P_{i} >= eps*eye(Nx),...
        SgmF_{i} >= eps*eye(Ny+Ny),...
		H_{i} <= [[2*(alp + 1e-6)*P_{i}, eps*ones(Nx,Nw)];[eps*ones(Nx,Nw)', eps*eye(Nw)]]...%(alp + 1e-6) to de-homogenize LMI
	    ];
end

%Objective function	
obj=alp;
end

function [vars,lmis,obj]=LMI_PI(Poly,gL2,KPlb,KPub,TIlb,TIub)
%Parameters ---------------------------------------------------------
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nw=Poly.Nw;		%Exogenous input variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces
A=Poly.A;		%Autonomous system matrix polytope
Bu=Poly.Bu;		%System maniputaled inputs polytope
Bw=Poly.Bw;		%System exogenous inputs polytope
Cy=Poly.Cy;		%Controlled output-state relationship polytope
Tl=Poly.Tl;     %Left-multiplier matrix of equivalence relation
Tr=Poly.Tr;     %Right-multiplier matrix of equivalence relation
R=Poly.R;      	%Weighting matrix of equivalence relation

%SDP variables ------------------------------------------------------
KP=sdpvar(Ny,Ny,'diagonal');
KI=sdpvar(Ny,Ny,'diagonal');

%Augmented propotional gain
if Ny<Nx
	XP=sdpvar(Ny,Nx-Ny,'full');
	YP=sdpvar(Nx-Ny,Ny,'full');
	ZP=sdpvar(Nx-Ny,Nx-Ny,'full');
	Kprp=[[KP,XP];[YP,ZP]];
else
	Kprp=KP;
end

%Slack variable for quadratic form disguise
Sgm_FP=sdpvar(Ny,Ny,'diagonal');
epsl_FP=sdpvar(1);

for i=1:Nl
%Lyapunov matrix proportional term
	Gprp{i}=(Tl{i}*Kprp*Tr{i}') + (Tl{i}*Kprp*Tr{i}')';
	Pprp{i}=Gprp{i}/2;

%Entire Lyapunov matrix
	P_{i}= [Pprp{i} Bu{i}
			Bu{i}'  eye(Ny)];
end

%Eigenvalue
alp=sdpvar(1);

%Storing SDP variables
vars.KP=KP; vars.KI=KI;
vars.Kprp=Kprp; vars.Pprp=Pprp;
vars.P_=P_;
vars.alp=alp;
vars.Sgm_FP=Sgm_FP;
vars.epsl_FP=epsl_FP;

%LMI system constraints ---------------------------------------------
lmis=[...
	KP >= KPlb, KP <= KPub,...
	KI >= TIub\KP, KI <= TIlb\KP,...
	Sgm_FP >= eps*eye(Ny),...
	epsl_FP <= 1 ...		
	];

for i=1:Nl
% Slack variable LMI constraint
	MP=KPub'*R{i}*KPub; sqrtMP=sqrtm(MP);
    pencilP=(sqrtm(R{i})*(KP-epsl_FP*KPub)/sqrtMP);
    SgmF_{i}=[[epsl_FP*eye(Ny)-Sgm_FP, pencilP'];
              [pencilP               , eye(Ny) ]];

%F_(KP,KI):
	FP=...	%FP with quadratic form disguise using slack variable Sgm_FP
		( (A{i}'*Gprp{i}) + (A{i}'*Gprp{i})' )/2 ...
		+ ( (Bu{i}*KI*Cy{i})+(Bu{i}*KI*Cy{i})' ) ...
        -2*(sqrtMP*Cy{i})'*Sgm_FP*(sqrtMP*Cy{i});

    FI=Bu{i}'*A{i} + KI*Cy{i} - 2*R{i}*KP*Cy{i};

%H_(KP,KI):    
    HP=FP + (Cy{i}'*Cy{i} + eps*ones(size(Cy{i}'*Cy{i},1))); %adding an almost null matrix to handle numerical issues
    SP=Bw{i}'*Gprp{i}/2;
    SI=Bw{i}'*Bu{i};

    H_{i}=[[HP, FI'    , SP'           ];...
           [FI, -2*R{i}, SI'           ];...
           [SP, SI     , -gL2^2*eye(Nw)]];

    lmis=[lmis,...
	    P_{i} >= eps*eye(Nx+Ny),...
		SgmF_{i} >= eps*eye(Ny+Ny),...
        H_{i} <= [[2*(alp + 1e-6)*P_{i}, eps*ones(Nx+Ny,Nw)];[eps*ones(Nx+Ny,Nw)', eps*eye(Nw)]]...%(alp + 1e-6) to de-homogenize LMI
	    ];
end

%Objective function	
obj=alp;
end

function [vars,lmis,obj]=LMI_PID(Poly,gL2,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub)
%Parameters ---------------------------------------------------------
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nw=Poly.Nw;		%Exogenous input variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces
A=Poly.A;		%Autonomous system matrix polytope
Bu=Poly.Bu;		%System maniputaled inputs polytope
Bw=Poly.Bw;		%System exogenous inputs polytope
Cy=Poly.Cy;		%Controlled output-state relationship polytope
Tl=Poly.Tl;     %Left-multiplier matrix of equivalence relation
Tr=Poly.Tr;     %Right-multiplier matrix of equivalence relation
R=Poly.R;      	%Weighting matrix of equivalence relation

%SDP variables ------------------------------------------------------
KP=sdpvar(Ny,Ny,'diagonal');
KI=sdpvar(Ny,Ny,'diagonal');
KD=sdpvar(Ny,Ny,'diagonal');

%Augmented propotional gain
if Ny<Nx
	XP=sdpvar(Ny,Nx-Ny,'full');
	YP=sdpvar(Nx-Ny,Ny,'full');
	ZP=sdpvar(Nx-Ny,Nx-Ny,'full');
	Kprp=[[KP,XP];[YP,ZP]];
else
	Kprp=KP;
end

%Slack variable for quadratic form disguise
MD=(TDub*KPub)'*PhiD*(TDub*KPub); sqrtMD=sqrtm(MD);
Sgm_Pprp=sdpvar(Ny,Ny,'diagonal'); 
epsl_Pprp=sdpvar(1); 

Sgm_FP=sdpvar(Ny,Ny,'diagonal');
epsl_FP=sdpvar(1);

for i=1:Nl
%Lyapunov matrix proportional term
	Gprp{i}=(Tl{i}*Kprp*Tr{i}') + (Tl{i}*Kprp*Tr{i}')';

	Pprp{i}=Gprp{i}/2 ...	%Pprp with quadratic form disguise using slack variable Sgm_Pprp
		- ( (Bu{i}*KD*PhiD*Cy{i}) + (Bu{i}*KD*PhiD*Cy{i})' )/2 ...
		+ (sqrtMD*Cy{i})'*Sgm_Pprp*(sqrtMD*Cy{i});

%Entire Lyapunov matrix
	P_{i}= [Pprp{i} 			Bu{i}		(Bu{i}'-KD*Cy{i})'	
			Bu{i}'  			eye(Ny)		zeros(Ny)
			(Bu{i}'-KD*Cy{i})	zeros(Ny)	eye(Ny)/PhiD];
end

%Eigenvalue
alp=sdpvar(1);

%Storing SDP variables
vars.KP=KP; vars.KI=KI; vars.KD=KD;
vars.Kprp=Kprp; vars.Pprp=Pprp;
vars.P_=P_;
vars.alp=alp;
vars.Sgm_Pprp=Sgm_Pprp; 
vars.Sgm_FP=Sgm_FP;
vars.epsl_Pprp=epsl_Pprp; 
vars.epsl_FP=epsl_FP;

%LMI system constraints ---------------------------------------------
% Slack variable LMI constraint
pencilD=(sqrtm(PhiD)*(KD-epsl_FP{i}*(KPub*TDub))/sqrtMD);
SgmPprp_=[[epsl_Pprp*eye(Ny)-Sgm_Pprp, pencilD'];
          [pencilD                   , eye(Ny) ]];

lmis=[...
	KP >= KPlb, KP <= KPub,...
	KI >= TIub\KP, KI <= TIlb\KP,...
	KD >= TDlb*KP, KD <= TDub*KP,...
    SgmPprp_ >= eps*eye(Ny+Ny),...
    Sgm_Pprp >= eps*eye(Ny),...
    epsl_Pprp <= 1,...	
	Sgm_FP >= eps*eye(Ny),...
	epsl_FP <= 1 ...	
	];

for i=1:Nl
%Slack variable LMI constraint
	MP=KPub'*R{i}*KPub; sqrtMP=sqrtm(MP);
    pencilP=(sqrtm(R{i})*(KP-epsl_FP*KPub)/sqrtMP);
    SgmF_{i}=[[epsl_FP*eye(Ny)-Sgm_FP, pencilP'];
              [pencilP               , eye(Ny) ]];

%F_(KP,KI,KD):
	FP=...	%FP with quadratic form disguise using slack variable Sgm_FP
		( (A{i}'*Gprp{i}) + (A{i}'*Gprp{i})' )/2 ...
		+ ( (Bu{i}*KD*PhiD*Cy{i}*A{i}) + (Bu{i}*KD*PhiD*Cy{i}*A{i})' )/2 ...
        - ( (A{i}'*Bu{i}*KD*PhiD*Cy{i}) + (A{i}'*Bu{i}*KD*PhiD*Cy{i})' )/2 ...
		+ ( (Bu{i}*KI*Cy{i})+(Bu{i}*KI*Cy{i})' ) ...
		-2*(sqrtMP*Cy{i})'*Sgm_FP*(sqrtMP*Cy{i});

    FI=Bu{i}'*A{i} + KI*Cy{i} - 2*R{i}*KP*Cy{i};

	FD=Bu{i}'*A{i} - PhiD*(Bu{i}' - KD*Cy{i}) - 2*R{i}*KP*Cy{i};

%H_(KP,KI,KD):
    HP=FP + (Cy{i}'*Cy{i} + eps*ones(size(Cy{i}'*Cy{i},1)));  %adding an almost null matrix to handle numerical issues

	SP=Bw{i}'*(Gprp{i}/2 - ( (Bu{i}*KD*PhiD*Cy{i}) - (Bu{i}*KD*PhiD*Cy{i})' )/2);
	SI=Bw{i}'*Bu{i};
    SD=Bw{i}'*Bu{i};

    H_{i}=[[HP, FI'    , FD'              , SP'           ];...
           [FI, -2*R{i}, -2*R{i}          , SI'           ];...
           [FD, -2*R{i}, -2*(R{i}+eye(Ny)), SD'           ];...
           [SP, SI     , SD               , -gL2^2*eye(Nw)]];

    lmis=[lmis,...
	    P_{i} >= eps*eye(Nx+Ny+Ny),...
		SgmF_{i} >= eps*eye(Ny+Ny),...
		H_{i} <= [[2*(alp + 1e-6)*P_{i}, eps*ones(Nx+Ny+Ny,Nw)];[eps*ones(Nx+Ny+Ny,Nw)', eps*eye(Nw)]]...%(alp + 1e-6) to de-homogenize LMI
	    ];  
end

%Objective function	
obj=alp;
end