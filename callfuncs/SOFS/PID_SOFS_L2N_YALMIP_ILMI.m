function [KP,KI,KD,TI,TD,alp,P_,iter]=...
	PID_SOFS_L2N_YALMIP_ILMI(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,law)
%	Semidefinite programming (SDP) problem to obtain feedback PID 
%	controller gain which maximizes system decay rate with upper
%	bound on Hinf norm by L2 function norms using YALMIP interface. 
%	SDP seeks to maximize system decay rate solving the generalized 
%	eigenvalue problem (GEVP) using YALMIP bisection built-in function.
%
%	The PID tuning is cast as a static output feedback (SOF) stabili-
%	zation problem with P_ and K_ as matrix variables. The bilinear
%	terms in P_ and K_ are handled by the iterative LMI approach of 
%	Wang et. al. (2008) PID Control for Multivariable Processes, 
%	Section 7.5.
%
% 	[KP,KI,KD,TI,TD,alp,P_]=...
%		PID_SOFS_L2N_YALMIP_ILMI(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,law)
%
%---- Outputs -------------------------------------------------------
%   KP: 	Proportional action gain matrix
%   KI: 	Integral action gain matrix
%   KD: 	Derivative action gain matrix
%   TI: 	Integral time matrix
%   TD: 	Derivative time matrix
%   alp: 	Decay rate
%   P_: 	Augmented system Lyapunov matrix
%	iter:	Number of LMI iterations
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
Nx=Poly.Nx; 	%State variables size
Nu=Poly.Nu; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Tk=sign(Tk); 	%Control action structure normalization

%Solve LMI problem by iterative LMI approach ------------------------
switch law
    case 'P'
A_=Poly.A0;
Bu_=Poly.Bu0*Tk;
    
	case 'PI'
A_ = [Poly.A0 				zeros(Nx,Ny)
	  Poly.Cy0  			zeros(Ny,Ny)];

Bu_= [Poly.Bu0*Tk
	  zeros(Ny,Nu)];

    case 'PID'
A_ = [Poly.A0               zeros(Nx,Ny)	zeros(Nx,Ny)
      Poly.Cy0              zeros(Ny)		zeros(Ny)
	  PhiD*Poly.Cy0*Poly.A0 zeros(Ny)		-PhiD*eye(Ny)];

Bu_= [Poly.Bu0*Tk
	 zeros(Ny,Nu)
	 PhiD*Poly.Cy0*Poly.Bu0*Tk];
	 
end
% [Remark]: The step 1 algorithm do not work if A_ have null eigenvalues. 
% Thus, it is needed to add an almost null negative  diagonal to A_ or 
% directly guess a P_ matrix, dropping guess Ricatti equation guess.
Q0_=eye(size(A_,1));
[P_,~,~]=care((A_-1e-8*eye(size(A_,1))),Bu_,Q0_);	%Step 1
% P_=1e-8*eye(size(A_,1));

%YALMIP SDP setting options
% optOP1=sdpsettings(...
%     'solver','bisection',...
%     'bisection.absgaptol',1e-5, ...
%     'bisection.relgaptol',1e-3,...
%     'bisection.solver','sedumi',...
% 	'sedumi.eps',1e-8,...
%     'sedumi.cg.qprec',1,'sedumi.cg.maxiter',49,...
%     'sedumi.stepdif',0,...'sedumi.stepdif',2,...'sedumi.stepdif',1,...
%     'sedumi.free',0,...'sedumi.free',1,...
%     'sedumi.sdp',1,...'sedumi.sdp',0,...
%     'verbose',1);
optOP1=sdpsettings(...
    'solver','bisection',...
    'bisection.absgaptol',1e-5, ...
    'bisection.relgaptol',1e-3,...
    'bisection.solver','sdpt3',...
    'verbose',1);
optOP2=sdpsettings(...
    'solver','sdpt3',...
    'verbose',0);	

iter=0; maxiter=100; alp=[];
while iter<maxiter
	iter=iter+1;
	alp0=alp;
	X_=value(P_);
	disp(strjoin({'----- Iteration #',num2str(iter),' --------------------------------------------------'},''));
	disp('----- solving gevp for SOFHinf --------------------------------------');
	[vars,lmis,obj]=LMIOP1(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,X_,law);   		%Step 2
	
    stat=mybisection(lmis,obj,optOP1,alp);
	
    if stat.problem==1
		warning(stat.info);
		break;
	end
	% K_=value(vars.K_);
    P_=value(vars.P_);

% [Remark]: The step 3 algorithm was modified as below to stop only
% when alp achieves convergence. The original algorithm stops in the first 
% alp<0 value found. However, this lead to non-optimal system decay rate tunings.
    alp=value(obj);
    res=abs((alp-alp0)/alp);
    disp(strjoin({'----- alp* = ',num2str(alp),' alp0* = ',num2str(alp0),' res = ',num2str(res)},''));
    if res<1e-3	%Step 3
	    break;
    end
%{
% [Remark]: Updating P_ by P_* lead to oscilating steps into iterative LMI 
% approach. After disabling this step, alp monotocally decreases after each 
% iteration step, reducing algorithm computational effort
    disp('----- Minimizing tr(P) for next iteration ----------------------------');
	[vars,lmis,obj]=LMIOP2(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,X_,alp,law);	%Step 4
    optOP2.verbose=0;
	stat=optimize(lmis,obj,optOP2);

% [Remark]: The step 4 SOFHinf algorithm was modified as below. 
% P_ is updated by P_* only if the min tr(P) problem is feasible.
% It was not clear which value of alp+dalp must be given to achieve
% feasibility of min tr(P). Issues about feasibility of min tr(P) problem
% is mentioned in the paper where originally the algorithm was developed 
% (Cao et al, 1998).
    if stat.problem==1
	    disp('----- P_ not updated by P_* because tr(P) problem is unfeasible -----');
    else
	    % K_=value(K_);
	    P_=value(P_);
	    trP=trace(P_);
	    disp(strjoin({'----- P_ updated by P_*. tr(P) = ',num2str(trP)},''));
    end
%}
    if norm(X_*Bu_-P_*Bu_)<1e-3	%Step 5
        warning('It cannot be decided by this algorithm whether SOFHinf problem is solvable');
        break;
    end
end
if iter==maxiter
    warning('Maximum number of LMI iterations achieved');
end
%--------------------------------------------------------------------
%Output solution ----------------------------------------------------
if stat.problem==1 || alp>eps
	KP=[];
	KI=[]; TI=[];
	KD=[]; TD=[];
	alp=[];
	P_=[];
else
    alp=value(vars.alp);
    KI=[]; TI=[];
    KD=[]; TD=[];	
	
    switch law
        case 'P'
    KP=value(vars.KP);

    %KP unpermutation by Tk
    KP=Tk*KP;
	
        case 'PI'
    KP=value(vars.KP);
    KI=value(vars.KI); TI=KP/KI;

    %KP, KI, TI unpermutation by Tk
    KP=Tk*KP;
    KI=Tk*KI; TI=abs(Tk)*TI;
	
        case 'PID'
    KP=value(vars.KP);
    KI=value(vars.KI); TI=KP/KI;
    KD=value(vars.KD); TD=KD/KP;	
    
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

function [vars,lmis,obj]=...
	LMIOP1(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,X_,law)

yalmip('clear');

%Parameters ---------------------------------------------------------
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nw=Poly.Nw;		%Exogenous input variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces

A_=cell(Nl,1); Bu_=cell(Nl,1); Bw_=cell(Nl,1); Cy_=cell(Nl,1);
for i=1:Nl
    switch law
        case 'P'
A_{i}=Poly.A{i};
Bu_=Poly.Bu{i}*Tk;
Bw_=Poly.Bw{i};
Cy_{i}=Poly.Cy{i};
    
	    case 'PI'
A_{i} = [Poly.A{i} 			zeros(Nx,Ny)
		 Poly.Cy{i}  		zeros(Ny,Ny)];

Bu_{i}= [Poly.Bu{i}*Tk
		 zeros(Ny)];

Bw_{i}= [Poly.Bw{i}
		 zeros(Ny,Nw)];
	 
Cy_{i} = [Poly.Cy{i}     	zeros(Ny,Ny)
		  zeros(Ny,Nx) 		eye(Ny,Ny)];	

        case 'PID'
A_{i} = [Poly.A{i}              	zeros(Nx,Ny)	zeros(Nx,Ny)
		 Poly.Cy{i}              	zeros(Ny)		zeros(Ny)
		 PhiD*Poly.Cy{i}*Poly.A{i} 	zeros(Ny)		-PhiD*eye(Ny)];

Bu_{i}= [Poly.Bu{i}*Tk
		 zeros(Ny)
		 PhiD*Poly.Cy{i}*Poly.Bu{i}*Tk];

Bw_{i}= [Poly.Bw{i}
		 zeros(Ny,Nw)
		 PhiD*Poly.Cy{i}*Poly.Bw{i}];	 
	 
Cy_{i}= [Poly.Cy{i}        	zeros(Ny)       zeros(Ny)
		 zeros(Ny,Nx)   	eye(Ny)			zeros(Ny)
		 zeros(Ny,Nx)		zeros(Ny)		eye(Ny)];
	 
    end
end
%SDP variables ------------------------------------------------------

% [Remark]: In the ILMI approach of Wang et. al. (2008), the eigenvalue term
% 2αP_ appears only in the first diagonal block matrix, for which B(z) is
% structurally semi-definite in the linear-fractional constraints A(z) <
% αB(z). This showed as a drawback for ILMI approach because SDP do not handle
% well strctural semi-definite matrices. Thus, a slack variable Sgm_ was created
% to handle this issue, with only the first diagonal block matrix considered as 
% linear-fractional constraint.

switch law
    case 'P'	%Lyapunov matrix for P control
KP=sdpvar(Ny,Ny,'diagonal');
K_=KP;
P_=sdpvar(Nx);
Sgm_=sdpvar(Nx);		 %slack variable

    case 'PI'	%Lyapunov matrix for PI control
KP=sdpvar(Ny,Ny,'diagonal');
KI=sdpvar(Ny,Ny,'diagonal');
K_=[KP,KI];
P_=sdpvar(Nx+Ny);
Sgm_=sdpvar(Nx+Ny);		 %slack variable

    case 'PID'	%Lyapunov matrix for PID control
KP=sdpvar(Ny,Ny,'diagonal');
KI=sdpvar(Ny,Ny,'diagonal');
KD=sdpvar(Ny,Ny,'diagonal');
K_=[KP,KI,KD];
P_=sdpvar(Nx+Ny+Ny);
Sgm_=sdpvar(Nx+Ny+Ny);   %slack variable
end

%Eigenvalue
alp=sdpvar(1);

%Storing SDP variables
vars.P_=P_;
vars.Sgm_=Sgm_;
vars.alp=alp;

switch law
    case 'P'	%Lyapunov matrix for P control
vars.KP=KP; vars.K_=K_;
    
    case 'PI'	%Lyapunov matrix for PI control
vars.KP=KP; vars.KI=KI; vars.K_=K_;

    case 'PID'	%Lyapunov matrix for PID control
vars.KP=KP; vars.KI=KI; vars.KD=KD; vars.K_=K_;
end

%LMI system constraints ---------------------------------------------
switch law
    case 'P'	%Lyapunov matrix for P control
lmis=[...
    P_ >= eps*eye(Nx),...
    Sgm_ <= eps*eye(Nx),...
	KP >= KPlb, KP <= KPub...
	];
    case 'PI'	%Lyapunov matrix for PI control
lmis=[...
    P_ >= eps*eye(Nx+Ny),...
    Sgm_ <= eps*eye(Nx+Ny),...
	KP >= KPlb, KP <= KPub,...
	KI >= TIub\KP, KI <= TIlb\KP...
	];
    case 'PID'	%Lyapunov matrix for PID control
lmis=[...
    P_ >= eps*eye(Nx+Ny+Ny),...
    Sgm_ <= eps*eye(Nx+Ny+Ny),...
	KP >= KPlb, KP <= KPub,...
	KI >= TIub\KP, KI <= TIlb\KP...
	KD >= TDlb*KP, KD <= TDub*KP...
	];
end

Nx_=size(P_,1);
for i=1:Nl
% H{i} < 0
% H{i}(P_,K_)
	% = A_'*P_ + P_*A_ + Cy'*Cy + P_*((Bw_/gL2^2)*Bw_')*P_ - P_*Bu_*K_*Cy_ - Cy_'*K_'*Bu_'*P_ < 2*alp*P_
	% = A_'*P_ + P_*A_ + Cy'*Cy + P_*((Bw_/gL2^2)*Bw_')*P_ - P_*(Bu_*Bu_')*P_ + (Bu_'*P_-K_*Cy_)'*(Bu_'*P_-K_*Cy_) - (K_*Cy_)'*(K_*Cy_) < 2*alp*P_
% dropping - (K_*Cy_)'*(K_*Cy_) yields to:
	% = A_'*P_ + P_*A_ + Cy'*Cy - Psi_ + P_*((Bw_/gL2^2)*Bw_')*P_ + (Bu_'*P_-K_*Cy_)'*(Bu_'*P_-K_*Cy_) < 2*alp*P_
% By schur complements
% [A_'*P_ + P_*A_ + Cy'*Cy - Psi_	P_*Bw_		(Bu_'*P_-K_*Cy_)'] < [2*alp*P_  0   0  ]
% [Bw_'*P_ 				            -gL2^2		0 				 ]   [0         eps 0  ]
% [(Bu_'*P_-K_*Cy_)		            0			-I 				 ]   [0         0	eps]
% Psi_ = X_'*(Bu_*Bu_')*P_ + P_*(Bu_*Bu_')*X_ - X_'*(Bu_*Bu_')*X_
    H11=...
	    ( A_{i}'*P_+P_*A_{i} ) ...
	    - ( (X_'*(Bu_{i}*Bu_{i}')*P_) + (X_'*(Bu_{i}*Bu_{i}')*P_)' ) + ( (Bu_{i}'*X_)'*(Bu_{i}'*X_) ) ...
        + Cy_{i}(1:Ny,:)'*Cy_{i}(1:Ny,:);
    
    H21=Bu_{i}'*P_ - K_*Cy_{i};
    
    H31=Bw_{i}'*P_;
    
    H_{i}=[[H11+Sgm_, H21'        , H31'          ];...
	       [H21     , -eye(Ny)    , zeros(Ny,Nw)  ];...
	       [H31     , zeros(Nw,Ny), -gL2^2*eye(Nw)]];
    
    lmis=[lmis,...
        H_{i} <= eps*eye(Nx_+Ny+Nw)...
        H11-Sgm_ <= 2*(alp+1e-6)*P_;
	    ];
end

%Objective function	
obj=alp;
end

function [vars,lmis,obj]=...
	LMIOP2(Poly,gL2,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,X_,alp,law)

yalmip('clear');

%Parameters ---------------------------------------------------------
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nw=Poly.Nw;		%Exogenous input variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces

A_=cell(Nl,1); Bu_=cell(Nl,1); Bw_=cell(Nl,1); Cy_=cell(Nl,1);
for i=1:Nl
    switch law
        case 'P'
A_{i}=Poly.A{i};
Bu_=Poly.Bu{i}*Tk;
Bw_=Poly.Bw{i};
Cy_{i}=Poly.Cy{i};
    
	    case 'PI'
A_{i} = [Poly.A{i} 			zeros(Nx,Ny)
		 Poly.Cy{i}  		zeros(Ny,Ny)];

Bu_{i}= [Poly.Bu{i}*Tk
		 zeros(Ny)];

Bw_{i}= [Poly.Bw{i}
		 zeros(Ny,Nw)];
	 
Cy_{i} = [Poly.Cy{i}     	zeros(Ny,Ny)
		  zeros(Ny,Nx) 		eye(Ny,Ny)];	

        case 'PID'
A_{i} = [Poly.A{i}              	zeros(Nx,Ny)	zeros(Nx,Ny)
		 Poly.Cy{i}              	zeros(Ny)		zeros(Ny)
		 PhiD*Poly.Cy{i}*Poly.A{i}  zeros(Ny)		-PhiD*eye(Ny)];

Bu_{i}= [Poly.Bu{i}*Tk
		 zeros(Ny)
		 PhiD*Poly.Cy{i}*Poly.Bu{i}*Tk];

Bw_{i}= [Poly.Bw{i}
		 zeros(Ny,Nw)
		 PhiD*Poly.Cy{i}*Poly.Bw{i}];	
	 
Cy_{i}= [Poly.Cy{i}        	zeros(Ny)       zeros(Ny)
		 zeros(Ny,Nx)   	eye(Ny)			zeros(Ny)
		 zeros(Ny,Nx)		zeros(Ny)		eye(Ny)];
	 
    end
end
%SDP variables ------------------------------------------------------
switch law
    case 'P'	%Lyapunov matrix for P control
KP=sdpvar(Ny,Ny,'diagonal');
K_=KP;
P_=sdpvar(Nx);
    
    case 'PI'	%Lyapunov matrix for PI control
KP=sdpvar(Ny,Ny,'diagonal');
KI=sdpvar(Ny,Ny,'diagonal');
K_=[KP,KI];
P_=sdpvar(Nx+Ny);

    case 'PID'	%Lyapunov matrix for PID control
KP=sdpvar(Ny,Ny,'diagonal');
KI=sdpvar(Ny,Ny,'diagonal');
KD=sdpvar(Ny,Ny,'diagonal');
K_=[KP,KI,KD];
P_=sdpvar(Nx+Ny+Ny);
end

%Storing SDP variables
switch law
    case 'P'	%Lyapunov matrix for P control
vars.KP=KP; vars.K_=K_;
vars.P_=P_; vars.alp=alp;
    
    case 'PI'	%Lyapunov matrix for PI control
vars.KP=KP; vars.KI=KI; vars.K_=K_;
vars.P_=P_; vars.alp=alp;

    case 'PID'	%Lyapunov matrix for PID control
vars.KP=KP; vars.KI=KI; vars.KD=KD; vars.K_=K_;
vars.P_=P_;
end

%LMI system constraints ---------------------------------------------
switch law
    case 'P'	%Lyapunov matrix for P control
lmis=[...
    P_ >= eps*eye(Nx),...
	KP >= KPlb, KP <= KPub...
	];
    case 'PI'	%Lyapunov matrix for PI control
lmis=[...
    P_ >= eps*eye(Nx+Ny),...
	KP >= KPlb, KP <= KPub,...
	KI >= TIub\KP, KI <= TIlb\KP...
	];
    case 'PID'	%Lyapunov matrix for PID control
lmis=[...
    P_ >= eps*eye(Nx+Ny+Ny),...
	KP >= KPlb, KP <= KPub,...
	KI >= TIub\KP, KI <= TIlb\KP...
	KD >= TDlb*KP, KD <= TDub*KP...
	];
end

Nx_=size(P_,1);
for i=1:Nl
% H{i} < 0
% H{i}(P_,K_)
	% = A_'*P_ + P_*A_ + Cy'*Cy + P_*((Bw_/gL2^2)*Bw_')*P_ - P_*Bu_*K_*Cy_ - Cy_'*K_'*Bu_'*P_ < 2*alp*P_
	% = A_'*P_ + P_*A_ + Cy'*Cy + P_*((Bw_/gL2^2)*Bw_')*P_ - P_*(Bu_*Bu_')*P_ + (Bu_'*P_-K_*Cy_)'*(Bu_'*P_-K_*Cy_) - (K_*Cy_)'*(K_*Cy_) < 2*alp*P_
% dropping - (K_*Cy_)'*(K_*Cy_) yields to:
	% = A_'*P_ + P_*A_ + Cy'*Cy - Psi_ -2*alp*P_ + P_*((Bw_/gL2^2)*Bw_')*P_ + (Bu_'*P_-K_*Cy_)'*(Bu_'*P_-K_*Cy_) < 0
% By schur complements
% [A_'*P_ + P_*A_ + Cy'*Cy - Psi_ - 2*alp*P_	P_*Bw_		(Bu_'*P_-K_*Cy_)'] < 0
% [Bw_'*P_ 										-gL2^2		0 				 ]
% [(Bu_'*P_-K_*Cy_)								0			-I 				 ]
% Psi_ = X_'*(Bu_*Bu_')*P_ + P_*(Bu_*Bu_')*X_ - X_'*(Bu_*Bu_')*X_
    H11=...
	    ( A_{i}'*P_+P_*A_{i} ) ...
	    - ( (X_'*(Bu_{i}*Bu_{i}')*P_) + (X_'*(Bu_{i}*Bu_{i}')*P_)' ) + ( (Bu_{i}'*X_)'*(Bu_{i}'*X_) ) ...
        + Cy_{i}(1:Ny,:)'*Cy_{i}(1:Ny,:) ...
        -2*alp*P_;
    
    H21=Bu_{i}'*P_ - K_*Cy_{i};
    
    H31=Bw_{i}'*P_;
    
    H_{i}=[[H11, H21'        , H31'          ];...
	       [H21, -eye(Ny)    , zeros(Ny,Nw)  ];...
	       [H31, zeros(Nw,Ny), -gL2^2*eye(Nw)]];
    
    lmis=[lmis,...
        H_{i} <= eps*eye(Nx_+Ny+Nw)...
	    ];
end

%Objective function	
obj=trace(P_);
end