close all;
T=200;
cGAMMA = 0.4;
cTAU   = 0.2;
cPHI   = 0.04;
cCHI   = 0.5;
cBETA  = 0.99;

SS=zeros(1,T);
II=zeros(1,T);
YY=zeros(1,T);
UU=zeros(1,T);
WW=zeros(1,T); 
RR=zeros(1,T);
DD=zeros(1,T);

% ============= 
cTAU_base   = cTAU;
cALPHA=zeros(1,T);

% Initial Condition
II(1) = 0.001;
SS(1) = 1 - II(1);
RR(1) = 0;
DD(1) = 0;
YY(1) = (SS(1)+RR(1))*(1-cALPHA(1)); 
UU(1) = YY(1) - cCHI*II(1);

%========= Dynamics
for t=2:T
	SS(t) = SS(t-1) - cGAMMA*(1-cALPHA(t-1))*SS(t-1)*II(t-1);	
	II(t) = II(t-1) + cGAMMA*(1-cALPHA(t-1))*SS(t-1)*II(t-1) - cTAU*II(t-1) - cPHI*II(t-1);	
	RR(t) = RR(t-1) + cTAU*II(t-1);	
	DD(t) = DD(t-1) + cPHI*II(t-1);	
	YY(t) = (SS(t)+RR(t))*(1-cALPHA(t)); 
	UU(t) = YY(t) - cCHI*DD(t);
end

%========= Steady state values
SSs = SS(T);
IIs = 0;
RRs = RR(T);
DDs = DD(T);
YYs = SSs+RRs;
UUs = YYs -cCHI*DDs;
WWs = UUs/(1-cBETA);

%========= Computing welfare
WW(T)=WWs;
for t=T-1:-1:1
	WW(t) = UU(t) +cBETA*WW(t+1);
end

figure(1)
subplot(2,3,1)
plot([1:T],SS,'k','LineWidth',1.5);
title('Susceptible','FontSize',16)
ylim([0 1])
hold on
subplot(2,3,2)
plot([1:T],II,'k','LineWidth',1.5);
title('Infectious','FontSize',16)
hold on
subplot(2,3,3)
plot([1:T],RR,'k','LineWidth',1.5);
title('Recovered','FontSize',16)
hold on
subplot(2,3,4)
plot([1:T],DD,'k','LineWidth',1.5);
title('Dead','FontSize',16)
ylim([0 1])
hold on
subplot(2,3,5)
plot([1:T],UU,'k','LineWidth',1.5);
title('Utility','FontSize',16)
ylim([0 1])
hold on
subplot(2,3,6)
plot([1:T],WW,'k','LineWidth',1.5);
title('Welfare','FontSize',16)
hold on
%legend('\tau=0.5','\tau=0.25','FontSize',16,'Location','Southeast')

figure(2)
plot(SS,II,'k','LineWidth',1.5);
xline((cTAU_base+cPHI)/cGAMMA,'k')
ylim([0 1.5*max(II)]) 
xlabel('S','FontSize',14)
ylabel('I','FontSize',14)
hold on
%legend('\tau=0.5','\tau=0.5','\tau=0.25','\tau=0.25','FontSize',16,'Location','Northeast')

set(figure(1),'PaperOrientation','Landscape');
set(figure(1),'PaperPosition',[0 0 2*11 2*8.5]);
%print(figure(1),'-dpng','sird_IRFs.png');

set(figure(2),'PaperOrientation','Landscape');
%print(figure(2),'-dpng','sird_phase.png');
