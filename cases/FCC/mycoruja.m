clear, close all
[R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,SPfcc,Temperature,Vapor,Comp1,Comp2,Comp3,Comp4,Comp5,Comp6,Comp7,Comp8,Comp9,Comp10,Liquid,Holdup,EnthalL,EnthalV,Ttrack,SPfrac,ProductsMB,LPG,MVfrac,ufra,xfra,xc,products,errord,xfil,Valvesfrac,time]=mydynamic;
% save('cases\FCC\results\data\FCC_Frac_Steady_lmi.mat');
% save('cases\FCC\results\data\FCC_Frac_Steady_ilmi.mat');
save('cases\FCC\results\data\FCC_Frac_Steady_nlp.mat'); %didn't worked. The numerical integration failed when using the tunings from nlp_FCC.mat 

xfcc=R5(end,2:end);
ufcc=[165,0,1,0,0,0,0.431260000000000,0,0,616,24.9000000000000,28,1250,98100,969];% LAST 6 COMPONENTS ARE THE SET POINTS FOR FCC. 
%Detailed definiton of each ufcc component can be found in FCC.m file
dist=[75,25,460.900000000000,0.900000000000000]; %DISTURBANCES
%Detailed definiton of each dist component can be found in FCC.m file
Flpg=LPG(end);
Tcondenser=Temperature(end,1);
MV=MVfrac(end,:)';
SP=SPfrac(end,:)';
Distillateini=xfra(20*3+1)/MV(1);
Xfilin=xfil(end,:)';
save('mydynamic_ss.mat','dist','Distillateini','errord','Flpg','MV','products','SP','Tcondenser','ufcc','ufra','xc','xfcc','Xfilin','xfra')
save('cases\FCC\results\data\mydynamic_ss.mat','dist','Distillateini','errord','Flpg','MV','products','SP','Tcondenser','ufcc','ufra','xc','xfcc','Xfilin','xfra');

clear, close all
[R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,SPfcc,Temperature,Vapor,Comp1,Comp2,Comp3,Comp4,Comp5,Comp6,Comp7,Comp8,Comp9,Comp10,Liquid,Holdup,EnthalL,EnthalV,Ttrack,SPfrac,ProductsMB,LPG,MVfrac,ufra,xfra,xc,products,errord,xfil,Valvesfrac,time]=mydynamic_SPtrack();
% save('cases\FCC\results\data\FCC_Frac_SPtrack_lmi.mat');
% save('cases\FCC\results\data\FCC_Frac_SPtrack_ilmi.mat');
save('cases\FCC\results\data\FCC_Frac_SPtrack_nlp.mat'); %didn't worked. The numerical integration failed when using the tunings from nlp_FCC.mat 

clear
[R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,SPfcc,Temperature,Vapor,Comp1,Comp2,Comp3,Comp4,Comp5,Comp6,Comp7,Comp8,Comp9,Comp10,Liquid,Holdup,EnthalL,EnthalV,Ttrack,SPfrac,ProductsMB,LPG,MVfrac,ufra,xfra,xc,products,errord,xfil,Valvesfrac,time]=mydynamic_DistRej();
% save('cases\FCC\results\data\FCC_Frac_DistRej_lmi.mat');
% save('cases\FCC\results\data\FCC_Frac_DistRej_ilmi.mat');
save('cases\FCC\results\data\FCC_Frac_DistRej_nlp.mat'); %didn't worked. The numerical integration failed when using the tunings from nlp_FCC.mat 

% myPlotall_SPtrack
% myPlotall_DistRej