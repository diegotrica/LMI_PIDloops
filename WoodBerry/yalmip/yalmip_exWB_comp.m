clear, close all

load('cases\WoodBerry\yalmip\results\data\yalmip_exWB_tb_rlt.mat','rlt');

tstart=0; tend=300; tstep=1; t=(tstart:tstep:tend)'; %time domain
yf=[0.1;0.1]; dw=0.245;
InputName={'w_{D,r} [wt. frac]','w_{B,r} [wt. frac]','f [lb/min]'};
OutputName={'w_D [wt. frac]','w_B [wt. frac]'};
Ny=size(yf,1); Nw=size(dw,1);

for Nrlt=1:3
    sys=ss(rlt(Nrlt).TW.PI.Gsp.A,[rlt(Nrlt).TW.PI.Gsp.B,rlt(Nrlt).TW.PI.Gload.B],rlt(Nrlt).TW.PI.Gsp.C,[]);
    yPI_TW=step(sys,t,stepDataOptions('StepAmplitude',[yf;dw]));
    
    sys=ss(rlt(Nrlt).ILMI.PI.Gsp.A,[rlt(Nrlt).ILMI.PI.Gsp.B,rlt(Nrlt).ILMI.PI.Gload.B],rlt(Nrlt).ILMI.PI.Gsp.C,[]);
    yPI_ILMI=step(sys,t,stepDataOptions('StepAmplitude',[yf;dw]));
    
    sys=ss(rlt(Nrlt).TW.PID.Gsp.A,[rlt(Nrlt).TW.PID.Gsp.B,rlt(Nrlt).TW.PID.Gload.B],rlt(Nrlt).TW.PID.Gsp.C,[]);
    yPID_TW=step(sys,t,stepDataOptions('StepAmplitude',[yf;dw]));
    
    sys=ss(rlt(Nrlt).ILMI.PID.Gsp.A,[rlt(Nrlt).ILMI.PID.Gsp.B,rlt(Nrlt).ILMI.PID.Gload.B],rlt(Nrlt).ILMI.PID.Gsp.C,[]);
    yPID_ILMI=step(sys,t,stepDataOptions('StepAmplitude',[yf;dw]));
    
    %Rendering figure
    figh=figure(Nrlt); pos=0;
    ytick{1}=(-1.5:0.5:1.5)'*yf(1);
    ytick{2}=(-0.4:0.8:4.4)'*yf(2);
    
    %Ploting results
    for i=1:Ny
	    for j=1:Ny+Nw
		    pos=pos+1;
		    ax=subplot(Ny,Ny+Nw,pos);
		    plot(ax,t,yPI_TW(:,i,j),'b-',t,yPI_ILMI(:,i,j),'b--',t,yPID_TW(:,i,j),'r-',t,yPID_ILMI(:,i,j),'r--');
            if i==1
			    title(ax,InputName{j},'FontName','Helvetica','FontSize',8);
		    end
		    xlabel(ax,'t [min]','FontName','Helvetica','FontSize',8,'HorizontalAlignment','center','Rotation',0);
		    set(ax,'XLim',[0 t(end)],'XTick',(0:60:t(end))','YLim',[ytick{i}(1) ytick{i}(end)],'YTick',ytick{i},'FontName','Helvetica','FontSize',8);
            hold on;
	    end
    end
    
    ax=subplot(Ny,Ny+Nw,1);
    ylabel(ax,OutputName{1},'FontName','Helvetica','FontSize',8);
    
    ax=subplot(Ny,Ny+Nw,1+Ny+Nw);
    ylabel(ax,OutputName{2},'FontName','Helvetica','FontSize',8);
    
    figh.Units='centimeters';
    figh.Color='white';
    figh.OuterPosition=[0 0 19 10.5];
    figh.Position=[0 0 19 10.5];
    
    %Set legend
    h=legend('PI-TW','PI-ILMI','PID-TW','PID-ILMI');
    h.Position=[0.40,0.955,0.2,0.05];
    h.Orientation='horizontal';
    h.FontName='Helvetica';
    h.FontSize=8;
    
    %Saving figure
    print(figh,rlt(Nrlt).data(1:end-4),'-dsvg');
    print(figh,rlt(Nrlt).data(1:end-4),'-depsc');
end