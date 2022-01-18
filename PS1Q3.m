close all;
T=200;
cGAMMA = 0.5;
cTAU   = 0.2;

SS=zeros(1,T);
II=zeros(1,T);
RR=zeros(1,T);

% Initial Condition
II(1) = 0.001;
SS(1) = 1 - II(1);
RR(1) = 0;

% Dynamics without vaccine
for t=2:T
    SS(t) = SS(t-1) - cGAMMA*SS(t-1)*II(t-1);	
	II(t) = II(t-1) + cGAMMA*SS(t-1)*II(t-1) - cTAU*II(t-1);	
	RR(t) = RR(t-1) + cTAU*II(t-1);	
end

X0    = cGAMMA/cTAU;
Sstar = 1/X0;
Rstar = 1-Sstar;

subplot(2,2,1)
plot([1:T],SS,'k','LineWidth',1.5);
hold on
plot([1:T],ones(1,T)*Sstar,'k','LineWidth',1);
title('Susceptible','FontSize',20)
subplot(2,2,2)
plot([1:T],II,'k','LineWidth',1.5);
hold on
title('Infectious','FontSize',20)
xlabel('Time')
subplot(2,2,3)
plot([1:T],RR,'k','LineWidth',1.5);
hold on
plot([1:T],ones(1,T)*Rstar,'k','LineWidth',1);
title('Recovered','FontSize',20)

SS_base=SS;
II_base=II;

% Introduce vaccine
V = zeros(1,T);
V(11:30) = 0.01;

% Dynamics with vaccine
for t=2:T
    SS(t) = SS(t-1) - cGAMMA*SS(t-1)*II(t-1) - V(t-1);	
	II(t) = II(t-1) + cGAMMA*SS(t-1)*II(t-1) - cTAU*II(t-1);	
	RR(t) = RR(t-1) + cTAU*II(t-1) + V(t-1);	
end

figure(1)
subplot(2,2,1)
plot([1:T],SS,'r','LineWidth',1.5);
hold on
title('Susceptible','FontSize',20)
subplot(2,2,2)
plot([1:T],II,'r','LineWidth',1.5);
hold on
title('Infectious','FontSize',20)
legend('No Vaccine','Vaccine','FontSize',16)
xlabel('Time')
subplot(2,2,3)
plot([1:T],RR,'r','LineWidth',1.5);
hold on
title('Recovered','FontSize',16)

figure(2)
plot(SS_base,II_base,'k','LineWidth',1.5);
hold on
plot(SS,II,'r','LineWidth',1.5);
xline(Sstar,'k','LineWidth',1)
ylim([0 1.5*max(max(II),max(II_base))]) 
xlabel('S','FontSize',14)
ylabel('I','FontSize',14)
legend('No Vaccine','Vaccine','FontSize',16)