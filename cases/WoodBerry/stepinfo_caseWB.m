function info=stepinfo_caseWB(KP,KI,KD,PhiD,sys,t,du,law)
Ny=size(KP,1); Nw=size(sys.B,2)-Ny;
switch law
	case 'PI'
Nx=size(sys.A,1)-Ny; 	
for i=1:Ny+Nw
Br_=[sys.B(1:Nx,1:Ny)*KP;
	 -eye(Ny)];
Bw_=[sys.B(1:Nx,1+Ny:end)
	 zeros(Ny,Nw)];
syscl=ss(sys.A-sys.B(:,1:Ny)*[KP,KI]*sys.C, [Br_,Bw_], sys.C(1:Ny,:), 0);
end

	case 'PID'
Nx=size(sys.A,1)-2*Ny;
for i=1:Ny+Nw
Br_=[sys.B(1:Nx,1:Ny)*KP;
	 -eye(Ny);
	 -PhiD*sys.C(1:Ny,1:Nx)*sys.B(1:Nx,1:Ny)*KP];
Bw_=[sys.B(1:Nx,1+Ny:end)
	 zeros(Ny,Nw)
	 -PhiD*sys.C(1:Ny,1:Nx)*sys.B(1:Nx,1+Ny:end)];
syscl=ss(sys.A-sys.B(:,1:Ny)*[KP,KI,KD]*sys.C, [Br_,Bw_], sys.C(1:Ny,:), 0);
end
end

y=step(syscl,t,stepDataOptions('StepAmplitude',du));
yf=diag(du);
for i=1:Ny
    for j=1:Ny+Nw
        info(i,j)=stepinfo(y(:,i,j),t,yf(i,j));
    end
end

dt=t(2:end)-t(1:end-1);
for i=1:Ny
    for j=1:Ny
        if i==j
			info(i,j).ITAE=sum(t(1:end-1).*abs(du(i)-y(1:end-1,i,j)).*dt);
        else
			info(i,j).ITAE=sum(t(1:end-1).*abs(y(1:end-1,i,j)).*dt);
        end
    end
    for j=1+Ny:Ny+Nw
		info(i,j).ITAE=sum(t(1:end-1).*abs(y(1:end-1,i,j)).*dt);
    end
end
end