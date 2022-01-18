close all;
T=21;

cH   = 1;
cD   = 1;
cA   = 2;
cws  = 1;
ca   = 0.4;

N   = zeros(1,T);
L   = zeros(1,T);
w   = zeros(1,T);
Y   = zeros(1,T);
YN  = zeros(1,T);


A   = 2*ones(1,T);
ETA = ones(1,T); 
D   = ones(1,T); 
PHI = ones(1,T); 
WS  = ones(1,T); 

% A transitory shock at time 1
%WS(1) = 0.5*cws; 

% A permanent shock at time 1
WS(1:T) = 0.5*cws; 

% Initial Steady state
cPHI = (1-ca)/cws*(cD/cH)^ca;

Nss = (cA*cPHI)^(1/ca);
Lss = cH*Nss;
wss = cws;
Yss = 2*cD^ca*Lss^(1-ca);
YNss = Yss/Nss;

% Initial Condition
N(1)  = Nss;                      
        
for t=1:T-1
	PHI(t)= (1-ca)/WS(t)*(D(t)/cH)^ca;
	N(t+1)= PHI(t)*A(t)*ETA(t)*N(t)^(1-ca);	
	L(t)  = cH*N(t);
	w(t)  = (1-ca)*A(t)*(D(t)/cH)^ca*N(t)^(-ca);
	Y(t)  = A(t)*D(t)^ca*L(t)^(1-ca);
	YN(t) = Y(t)/N(t);
end
t=T;
PHI(t)= (1-ca)/WS(t)*(D(t)/cH)^ca;
L(t)  = cH*N(t);
w(t)  = (1-ca)*A(t)*(D(t)/cH)^ca*N(t)^(-ca);
Y(t)  = A(t)*D(t)^ca*L(t)^(1-ca);
YN(t) = Y(t)/N(t);

figure(2)
subplot(2,2,1)
plot([0:T-1],N(1:T),'k','LineWidth',1.5);
hold on
plot([0:T-1],Nss*ones(1,T),'k');
xlim([0 T-1])
xlabel('Time')
title('Population','FontSize',14)
subplot(2,2,2)
plot([0:T-1],w(1:T),'k','LineWidth',1.5);
hold on
plot([0:T-1],wss*ones(1,T),'k');
xlim([0 T-1])
xlabel('Time')
title('wage','FontSize',14)
subplot(2,2,3)
plot([0:T-1],Y(1:T),'k','LineWidth',1.5);
hold on
plot([0:T-1],Yss*ones(1,T),'k');
title('output','FontSize',14)
xlim([0 T-1])
xlabel('Time')
subplot(2,2,4)
plot([0:T-1],WS(1:T),'k','LineWidth',1.5);
xlim([0 T-1])
title('Subsistence Wage','FontSize',14)
xlabel('Time')


