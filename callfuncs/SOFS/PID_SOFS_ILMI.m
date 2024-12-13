function [KP,KI,KD,TI,TD,alp,P_,iter]=...
	PID_SOFS_ILMI(Poly,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,law)
%	Semidefinite programming (SDP) problem to obtain feedback PID 
%	controller gain which maximizes system decay rate. SDP seeks to 
%	maximize system decay rate solving the generalized eigenvalue 
%	problem (GEVP) using the gevp command.
%
%	The PID tuning is cast as a static output feedback (SOF) stabili-
%	zation problem with P_ and K_ as matrix variables. The bilinear
%	terms in P_ and K_ are handled by the iterative LMI approach of 
%	Wang et. al. (2008) PID Control for Multivariable Processes, 
%	Section 7.5.
%
% 	[KP,KI,KD,TI,TD,alp,P_,iter]=...
%		PID_SOFS_ILMI(Poly,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,law)
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

%	Developed for t0he paper:
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
Nu=Ny;			%Input variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces
Tk=sign(Tk); 	%Control action structure normalization

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
%Solve SDP problem by iterative LMI approach ------------------------
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
% Thus, it is needed to add an almost null negative diagonal to A_ or 
% directly guess a P_ matrix, dropping guess Ricatti equation guess.
% Q0_=eye(size(A_,1));
% [P_,~,~]=care((A_-1e-8*eye(size(A_,1))),Bu_,Q0_);	%Step 1
P_=1e-8*eye(size(A_,1));

iter=0; maxiter=100; alp=[]; zOP1=[];
LMIOPT=[rtol,500,-1,0,0]; %LMI solver options
while iter<maxiter
	iter=iter+1;
	alp0=alp;
	X_=P_;
	disp(strjoin({'----- Iteration #',num2str(iter),' --------------------------------------------------'},''));
	disp('----- solving gevp for SOFS -----------------------------------------');
	LMISYS=LMIOP1(Poly,sK_,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,X_,law);  		%Step 2
	[alp,zOP1]=gevp(LMISYS,Nl,LMIOPT,alp,zOP1,-Inf);
	if isempty(alp)
		warning('Unfeasible problem. Check or modify constraints');
		break;
	end
	K_=dec2mat(LMISYS,zOP1,1);
	P_=dec2mat(LMISYS,zOP1,2);

% [Remark]: The step 3 algorithm was modified as below to stop only
% when alp achieves convergence. The original algorithm stops in the first 
% alp<0 value found. However, this lead to non-optimal system decay rate tunings.
	res=abs((alp-alp0)/alp);
	disp(strjoin({'----- alp* = ',num2str(alp),' alp0* = ',num2str(alp0),' res = ',num2str(res)},''));
	if res<1e-5	%Step 3
		break;
	end
%{
% [Remark]: Updating P_ by P_* lead to oscilating steps into iterative LMI 
% approach. After disabling this step, alp monotocally decreases after each 
% iteration step, reducing algorithm computational effort
	disp('----- Minimizing tr(P) for next iteration ---------------------------');
	LMISYS=LMIOP2(Poly,sK_,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,X_,alp,law);	%Step 4
	c=zeros(decnbr(LMISYS),1); idtrP_=diag(decinfo(LMISYS,2)); c(idtrP_)=1;
	zOP2=zOP1(1:size(K_(K_~=0),1)+size(P_(tril(P_)~=0),1));
	[trP,zOP2]=mincx(LMISYS,c,[LMIOPT(1:end-1),1],zOP2,0);

% [Remark]: The step 4 SOFHinf algorithm was modified as below. 
% P_ is updated by P_* only if the min tr(P) problem is feasible.
% It was not clear which value of alp+dalp must be given to achieve
% feasibility of min tr(P). Issues about feasibility of min tr(P) problem
% is mentioned in the paper where originally the algorithm was developed 
% (Cao et al, 1998).
	if ~isempty(zOP2)
		K_=dec2mat(LMISYS,zOP2,1);
		P_=dec2mat(LMISYS,zOP2,2);
		zOP1(1:size(K_(K_~=0),1)+size(P_(tril(P_)~=0),1))=zOP2;
		disp(strjoin({'----- P_ updated by P_*. tr(P) = ',num2str(trP)},''));
	else
		disp('----- P_ not updated by P_* because tr(P) problem is unfeasible -----');
	end
%}
	if norm(X_*Bu_-P_*Bu_)<1e-3	%Step 5
		warning('It cannot be decided by this algorithm whether SOF problem is solvable');
		break;
	end
end
if iter==maxiter
    warning('Maximum number of LMI iterations achieved');
end
%--------------------------------------------------------------------
%Output solution ----------------------------------------------------
if isempty(alp)
	KP=[];
	KI=[]; TI=[];
	KD=[]; TD=[];
	alp=[];
	P_=[];
else
	KI=[]; TI=[];
	KD=[]; TD=[];
	switch law
		case 'P'
	KP=K_(:,1:Ny);

	%KP unpermutation by Tk
	KP=Tk*KP;

		case 'PI'
	KP=K_(:,1:Ny);
	KI=K_(:,1+Ny:2*Ny);
	TI=KP/KI;

	%KP, KI, TI unpermutation by Tk
	KP=Tk*KP;
	KI=Tk*KI; TI=abs(Tk)*TI;

		case 'PID'
	KP=K_(:,1:Ny);
	KI=K_(:,1+Ny:2*Ny);
	TI=KP/KI;
	KD=K_(:,1+2*Ny:3*Ny);
	TD=KD/KP;

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

function LMISYS=...
	LMIOP1(Poly,sK_,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,X_,law)
%Parameters
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces

A_=cell(Nl,1); Bu_=cell(Nl,1); Cy_=cell(Nl,1);
for i=1:Nl
	switch law
		case 'P'
	A_{i}=Poly.A{i};
	Bu_=Poly.Bu{i}*Tk;
	Cy_{i}=Poly.Cy{i};
		
		case 'PI'
	A_{i} = [Poly.A{i} 			zeros(Nx,Ny)
			 Poly.Cy{i}  		zeros(Ny,Ny)];

	Bu_{i}= [Poly.Bu{i}*Tk
			 zeros(Ny)];
		 
	Cy_{i} = [Poly.Cy{i}     	zeros(Ny,Ny)
			  zeros(Ny,Nx) 		eye(Ny,Ny)];	

		case 'PID'
	A_{i} = [Poly.A{i}                  zeros(Nx,Ny)	zeros(Nx,Ny)
			 Poly.Cy{i}              	zeros(Ny)		zeros(Ny)
			 PhiD*Poly.Cy{i}*Poly.A{i}  zeros(Ny)		-PhiD*eye(Ny)];

	Bu_{i}= [Poly.Bu{i}*Tk
			 zeros(Ny)
			 PhiD*Poly.Cy{i}*Poly.Bu{i}*Tk]; 

	Cy_{i}= [Poly.Cy{i}        	zeros(Ny)       zeros(Ny)
			 zeros(Ny,Nx)   	eye(Ny)			zeros(Ny)
			 zeros(Ny,Nx)		zeros(Ny)		eye(Ny)];
		 
	end
end

%Initialize LMI system setting
setlmis([]);

% [Remark]: In the ILMI approach of Wang et. al. (2008), the eigenvalue term
% 2αP_ appears only in the first diagonal block matrix, for which B(z) is
% structurally semi-definite in the linear-fractional constraints A(z) <
% αB(z). This is a drawback for the MATLAB LMI Toolbox because B(z) must be
% positive-definite or forced by LMIs to be positive-definite. Thus, a slack
% variable Sgm_ was created to handle this issue, with only the first
% diagonal block matrix considered as linear-fractional constraint.

%LMI system variables
K_=lmivar(3,sK_);	%Augmented control gain

switch law
    case 'P'	%Lyapunov matrix for P control
P_=lmivar(1,[Nx 1]);
Sgm_=lmivar(1,[Nx 1]);		 %slack variable
    
    case 'PI'	%Lyapunov matrix for PI control
P_=lmivar(1,[Nx+Ny 1]);
Sgm_=lmivar(1,[Nx+Ny 1]);	 %slack variable

    case 'PID'	%Lyapunov matrix for PID control
P_=lmivar(1,[Nx+Ny+Ny 1]);
Sgm_=lmivar(1,[Nx+Ny+Ny 1]); %slack variable
end

KP=lmivar(3,sK_(:,1:Ny)); %Proportional gain

%LMI system constraints
neq=0;

%P_ > 0
neq=neq+1;
lmiterm([-neq 1 1 P_],1,1);

%Sgm_ < 0
neq=neq+1;
lmiterm([neq 1 1 Sgm_],1,1);

%|KP| > |KPlb|
neq=neq+1;
lmiterm([-neq 1 1 KP],1,1/2,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([neq 1 1 0],(KPlb+KPlb')/2); %Forcing symmetric form to not disp warning in command window

%|KP| < |KPub|
neq=neq+1;
lmiterm([neq 1 1 KP],1,1/2,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([-neq 1 1 0],(KPub+KPub')/2); %Forcing symmetric form to not disp warning in command window

switch law
	case 'PI'
KI=lmivar(3,sK_(:,1+Ny:2*Ny));		%Integral gain

%|KI| > |TIub\KP|
neq=neq+1;
lmiterm([-neq 1 1 KI],1/2,1,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([neq 1 1 KP],eye(Ny)/TIub,1/2,'s'); %Forcing symmetric form to not disp warning in command window
	
%|KI| < |TIlb\KP|
neq=neq+1;
lmiterm([neq 1 1 KI],1,1/2,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([-neq 1 1 KP],eye(Ny)/TIlb,1/2,'s'); %Forcing symmetric form to not disp warning in command window

    case 'PID'
KI=lmivar(3,sK_(:,1+Ny:2*Ny));		%Integral gain
KD=lmivar(3,sK_(:,1+2*Ny:3*Ny));	%Derivative gain

%|KI| > |TIub\KP|
neq=neq+1;
lmiterm([-neq 1 1 KI],1/2,1,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([neq 1 1 KP],eye(Ny)/TIub,1/2,'s'); %Forcing symmetric form to not disp warning in command window
	
%|KI| < |TIlb\KP|
neq=neq+1;
lmiterm([neq 1 1 KI],1,1/2,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([-neq 1 1 KP],eye(Ny)/TIlb,1/2,'s'); %Forcing symmetric form to not disp warning in command window

%|KD| > |TDlb*KP|
neq=neq+1;
lmiterm([-neq 1 1 KD],1,1/2,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([neq 1 1 KP],TDlb,1/2,'s'); %Forcing symmetric form to not disp warning in command window
	
%|KD| < |TDub*KP|
neq=neq+1;
lmiterm([neq 1 1 KD],1,1/2,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([-neq 1 1 KP],TDub,1/2,'s'); %Forcing symmetric form to not disp warning in command window
end

for i=1:Nl
% F{i} < 0
% F{i}(P_,K_)
	% = A_'*P_ + P_*A_ - P_*Bu_*K_*Cy_ - Cy_'*K_'*Bu_'*P_ < 2*alp*P_
	% = A_'*P_ + P_*A_ - P_*(Bu_*Bu_')*P_ + (Bu_'*P_-K_*Cy_)'*(Bu_'*P_-K_*Cy_) - (K_*Cy_)'*(K_*Cy_) < 2*alp*P_
% dropping - (K_*Cy_)'*(K_*Cy_) yields to:
	% = A_'*P_ + P_*A_ - Psi_ + (Bu_'*P_-K_*Cy_)'*(Bu_'*P_-K_*Cy_) < 2*alp*P_
% By schur complements
% [A_'*P_ + P_*A_ - Psi_	(Bu_'*P_-K_*Cy_)'] < [2*alp*P_	0  ]
% [(Bu_'*P_-K_*Cy_)		-I 				 ]   [0        	eps]
% Psi_ = X_'*(Bu_*Bu_')*P_ + P_*(Bu_*Bu_')*X_ - X_'*(Bu_*Bu_')*X_
	neq=neq+1;
	lmiterm([neq 1 1 Sgm_],1,1);
	lmiterm([neq 1 1 P_],A_{i}',1,'s');
	lmiterm([neq 1 1 P_],-X_'*(Bu_{i}*Bu_{i}'),1,'s');
	lmiterm([neq 1 1 0],(Bu_{i}'*X_)'*(Bu_{i}'*X_));
	lmiterm([neq 2 1 P_],Bu_{i}',1);
	lmiterm([neq 2 1 K_],1,-Cy_{i});
	lmiterm([neq 2 2 0],-1);

end

for i=1:Nl
%(LFC) A_'*P_ + P_*A_ - Psi_ - Sgm_ < 2*alp*P_	
	neq=neq+1;
	lmiterm([neq 1 1 P_],A_{i}',1,'s');
	lmiterm([neq 1 1 P_],-X_'*(Bu_{i}*Bu_{i}'),1,'s');
	lmiterm([neq 1 1 0],(Bu_{i}'*X_)'*(Bu_{i}'*X_));
	lmiterm([neq 1 1 Sgm_],-1,1);
	lmiterm([-neq 1 1 P_],1,1,'s');
end

%End and store LMI system
LMISYS=getlmis;
end

function LMISYS=...
	LMIOP2(Poly,sK_,Tk,PhiD,KPlb,KPub,TIlb,TIub,TDlb,TDub,X_,alp,law)
%Parameters
Nx=Poly.Nx; 	%State variables size
Ny=Poly.Ny;		%Output variables size
Nl=Poly.Nl;		%Number of polytope LTI system verteces

A_=cell(Nl,1); Bu_=cell(Nl,1); Cy_=cell(Nl,1);
for i=1:Nl
	switch law
		case 'P'
	A_{i}=Poly.A{i};
	Bu_=Poly.Bu{i}*Tk;
	Cy_{i}=Poly.Cy{i};
		
		case 'PI'
	A_{i} = [Poly.A{i} 			zeros(Nx,Ny)
			 Poly.Cy{i}  		zeros(Ny,Ny)];

	Bu_{i}= [Poly.Bu{i}*Tk
			 zeros(Ny)];
		 
	Cy_{i} = [Poly.Cy{i}     	zeros(Ny,Ny)
			  zeros(Ny,Nx) 		eye(Ny,Ny)];	

		case 'PID'
	A_{i} = [Poly.A{i}              	zeros(Nx,Ny)	zeros(Nx,Ny)
			 Poly.Cy{i}              	zeros(Ny)		zeros(Ny)
			 PhiD*Poly.Cy{i}*Poly.A{i} 	zeros(Ny)		-PhiD*eye(Ny)];

	Bu_{i}= [Poly.Bu{i}*Tk
			 zeros(Ny)
			 PhiD*Poly.Cy{i}*Poly.Bu{i}*Tk]; 

	Cy_{i}= [Poly.Cy{i}        	zeros(Ny)       zeros(Ny)
			 zeros(Ny,Nx)   	eye(Ny)			zeros(Ny)
			 zeros(Ny,Nx)		zeros(Ny)		eye(Ny)];
		 
	end
end

%Initialize LMI system setting
setlmis([]);
	
%LMI system variables
K_=lmivar(3,sK_);  					%Augmented control gain

switch law
    case 'P'                        %Lyapunov matrix for P control
P_=lmivar(1,[Nx 1]);
    
    case 'PI'                       %Lyapunov matrix for PI control
P_=lmivar(1,[Nx+Ny 1]);

    case 'PID'                      %Lyapunov matrix for PID control
P_=lmivar(1,[Nx+Ny+Ny 1]);
end

KP=lmivar(3,sK_(:,1:Ny));			%Proportional gain

%LMI system constraints
neq=0;

%P_ > 0
neq=neq+1;
lmiterm([-neq 1 1 P_],1,1);

%|KP| > |KPlb|
neq=neq+1;
lmiterm([-neq 1 1 KP],1,1/2,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([neq 1 1 0],(KPlb+KPlb')/2); %Forcing symmetric form to not disp warning in command window

%|KP| < |KPub|
neq=neq+1;
lmiterm([neq 1 1 KP],1,1/2,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([-neq 1 1 0],(KPub+KPub')/2); %Forcing symmetric form to not disp warning in command window

switch law
	case 'PI'
KI=lmivar(3,sK_(:,1+Ny:2*Ny));		%Integral gain

%|KI| > |TIub\KP|
neq=neq+1;
lmiterm([-neq 1 1 KI],1/2,1,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([neq 1 1 KP],eye(Ny)/TIub,1/2,'s'); %Forcing symmetric form to not disp warning in command window
	
%|KI| < |TIlb\KP|
neq=neq+1;
lmiterm([neq 1 1 KI],1,1/2,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([-neq 1 1 KP],eye(Ny)/TIlb,1/2,'s'); %Forcing symmetric form to not disp warning in command window

    case 'PID'
KI=lmivar(3,sK_(:,1+Ny:2*Ny));		%Integral gain
KD=lmivar(3,sK_(:,1+2*Ny:3*Ny));	%Derivative gain

%|KI| > |TIub\KP|
neq=neq+1;
lmiterm([-neq 1 1 KI],1/2,1,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([neq 1 1 KP],eye(Ny)/TIub,1/2,'s'); %Forcing symmetric form to not disp warning in command window
	
%|KI| < |TIlb\KP|
neq=neq+1;
lmiterm([neq 1 1 KI],1,1/2,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([-neq 1 1 KP],eye(Ny)/TIlb,1/2,'s'); %Forcing symmetric form to not disp warning in command window

%|KD| > |TDlb*KP|
neq=neq+1;
lmiterm([-neq 1 1 KD],1,1/2,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([neq 1 1 KP],TDlb,1/2,'s'); %Forcing symmetric form to not disp warning in command window
	
%|KD| < |TDub*KP|
neq=neq+1;
lmiterm([neq 1 1 KD],1,1/2,'s'); %Forcing symmetric form to not disp warning in command window
lmiterm([-neq 1 1 KP],TDub,1/2,'s'); %Forcing symmetric form to not disp warning in command window
end

for i=1:Nl
% F{i} < 0
% F{i}(P_,K_)
	% = A_'*P_ + P_*A_ - P_*Bu_*K_*Cy_ - Cy_'*K_'*Bu_'*P_ < 2*alp*P_
	% = A_'*P_ + P_*A_ - P_*(Bu_*Bu_')*P_ + (Bu_'*P_-K_*Cy_)'*(Bu_'*P_-K_*Cy_) - (K_*Cy_)'*(K_*Cy_) < 2*alp*P_
% dropping - (K_*Cy_)'*(K_*Cy_) yields to:
	% = A_'*P_ + P_*A_ - Psi_ -2*alp*P_ + (Bu_'*P_-K_*Cy_)'*(Bu_'*P_-K_*Cy_) < 0
% By schur complements
% [A_'*P_ + P_*A_ - Psi_ - 2*alp*P_	(Bu_'*P_-K_*Cy_)'] < 0
% [(Bu_'*P_-K_*Cy_)					-I 				 ]
% Psi_ = X_'*(Bu_*Bu_')*P_ + P_*(Bu_*Bu_')*X_ - X_'*(Bu_*Bu_')*X_
	neq=neq+1;
	lmiterm([neq 1 1 P_],A_{i}',1,'s');
	lmiterm([neq 1 1 P_],-alp,1,'s');
	lmiterm([neq 1 1 P_],-X_'*(Bu_{i}*Bu_{i}'),1,'s');
	lmiterm([neq 1 1 0],((X_'*Bu_{i})*(Bu_{i}'*X_)));
	lmiterm([neq 2 1 P_],Bu_{i}',1);
	lmiterm([neq 2 1 K_],1,-Cy_{i});
	lmiterm([neq 2 2 0],-1);
end

%End and store LMI system
LMISYS=getlmis;
end