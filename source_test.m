%% This code verifies that the correction terms W_1 and W_2 solve
% the corresponding PDE's. See stochasticVol.pdf for more info

%clear
clear all
close all
clc

%Parameters
T = 1 ;
rho = -0.10 ;
K = 30 ;
n = 2 ;


%functions
d2 = @(F,alpha,tau,strike) ( log( (1.0/strike)*F ) - 0.5*tau*(alpha.^2) )./(sqrt(tau)*alpha) ;
d1 = @(F,alpha,tau,strike) ( log( (1.0/strike)*F ) + 0.5*tau*(alpha.^2) )./(sqrt(tau)*alpha) ;
y = @(F,alpha,tau,strike) ( log( (1.0/strike)*F ) - 0.5*tau*(alpha.^2) )./(alpha);

N = @(x) normcdf(x) ;
N1 = @(x) normpdf(x) ;
N2 = @(x) -x.*normpdf(x) ;
N3 = @(x) (x.^2 - 1.0 ).*normpdf(x) ;
N4 = @(x) (-x.^3+3*x).*normpdf(x) ;
N5 = @(x) (x.^4 - 6*x.^2 + 3).*normpdf(x) ;
N6 = @(x) (-x.^5 + 10*x.^3 - 15*x).*normpdf(x) ;



W0 = @(F,a,t) F.*N(d1(F,a,T-t,K)) - K*N(d2(F,a,T-t,K)) 
W1 = @(F,a,t) 0.5*rho*K*a.*((T-t).*N2(d2(F,a,T-t,K)))
W2 = @(F,a,t) -(1/2)*K*a.*((1/3)*(a).*((T-t).^2).*N2(d2(F,a,T-t,K))...
    + (1/3)*((T-t).^(3/2)).*d2(F,a,T-t,K).*N2(d2(F,a,T-t,K))...
    -(1/6)*((T-t).^(3/2)).*N1(d2(F,a,T-t,K)))...
    -(1/2)*rho*rho*K*a.*((1/4)*(a).*((T-t).^2).*N4(d2(F,a,T-t,K))...
    + (1/4)*((T-t).^(3/2)).*d2(F,a,T-t,K).*N4(d2(F,a,T-t,K))...
    -(1/4)*((T-t).^(3/2)).*N3(d2(F,a,T-t,K)))
V = @(F,a,t) ((1/((n+2)*(n+3)))*(T-t) + (1/(n+3))*(y(F,a,T-t,K) + a.*(T-t)).^2).*((T-t).^(n+1-2)).*N4(d2(F,a,T-t,K))...
    -2*3*((1/((n+2)*(n+3)))*(y(F,a,T-t,K) + a.*(T-t))).*((T-t).^(n+1-3/2)).*N3(d2(F,a,T-t,K))...
    +2*3*2*(1/((n+1)*(n+2)*(n+3))).*((T-t).^(n)).*N2(d2(F,a,T-t,K))
V1 = @(F,t) -(1/(n+1))*((T-t).^(n-1)).*N4(d2(F,a,T-t,K))

%Plot
F = 20:1:60 ;
a = 0.20:.01:0.60;
t = 0 ;
[Fg, ag] = meshgrid(F,a) ;
z = W2(Fg,ag,t) ;
surf(Fg,ag,z)
title('W^2')

% Check source term
h = 0.00001;
dt = 0.000000001;
W1_AA = (1/(h*h))*(W1(Fg, ag+h, t) - 2*W1(Fg, ag, t)+ W1(Fg, ag-h, t)) ;
W2_FA = (1.0/(4*h*h))*(W2(Fg+h, ag+h, t) - W2(Fg+h, ag-h, t) - W2(Fg-h, ag+h, t)+W2(Fg-h, ag-h, t));
W0_AA = (1/(h*h))*(W0(Fg, ag+h, t) - 2*W0(Fg, ag, t)+ W0(Fg, ag-h, t)) ;
W1_FA = (1.0/(4*h*h))*(W1(Fg+h, ag+h, t) - W1(Fg+h, ag-h, t) - W1(Fg-h, ag+h, t)+W1(Fg-h, ag-h, t));
W0_FA = (1.0/(4*h*h))*(W0(Fg+h, ag+h, t) - W0(Fg+h, ag-h, t) - W0(Fg-h, ag+h, t)+W0(Fg-h, ag-h, t));
W2_A = (1.0/h)*(W2(Fg,ag+h,t)-W2(Fg,ag,t));
% V_FF = (1/(h*h))*(V(Fg+h,tg)-2*V(Fg,tg)+V(Fg-h,tg)) ;
% V_t = (1.0/dt)*(V(Fg,tg+dt)-V(Fg,tg));
% W2_FF = (1/(h*h))*(W2(Fg+h,tg)-2*W2(Fg,tg)+W2(Fg-h,tg)) ;
% W2_t = (1.0/dt)*(W2(Fg,tg+dt)-W2(Fg,tg));
% W1_FF = (1/(h*h))*(W1(Fg+h,tg)-2*W1(Fg,tg)+W1(Fg-h,tg)) ;
% W1_t = (1.0/dt)*(W1(Fg,tg+dt)-W1(Fg,tg));
% W0_FF = (1/(h*h))*(W0(Fg+h,tg)-2*W0(Fg,tg)+W0(Fg-h,tg)) ;
% W0_t = (1.0/dt)*(W0(Fg,tg+dt)-W0(Fg,tg));
% sourceV1 = ((T-tg).^(n-2)).*N4(d2(Fg,a,T-tg,K)) 
% source = ((T-tg).^(n-2)).*((a*(T-tg)+ y(Fg,a,T-tg,K)).^2).*N4(d2(Fg,a,T-tg,K));
source1 = -rho*K*ag.*N2(d2(Fg, ag, T-t, K)); 
source2 = (1/2)*K*ag.*((T-t)*ag+sqrt(T-t)*d2(Fg,ag,T-t,K)).*N2(d2(Fg,ag,T-t,K))...
    +(1/2)*rho*rho*K*ag.*((T-t)*ag+sqrt(T-t).*d2(Fg,ag,T-t,K)).*N4(d2(Fg,ag,T-t,K));
source_FA = -(1/4)*rho*K*(T-t)^(1/2)*ag.*((1/3)*(T-t)*ag+sqrt(T-t).*d2(Fg,ag,T-t,K)).*N3(d2(Fg,ag,T-t,K))...
    -(1/6)*rho*K*ag.*((y(Fg,ag,T-t,K)+(T-t)*ag).^2).*N4(d2(Fg,ag,T-t,K))...
    -(1/8)*(rho^3) * K *(T-t)^(1/2)*ag.*y(Fg,ag,T-t,K).*N5(d2(Fg,ag,T-t,K))...
    -(1/8)*(rho^3)*K*ag.*((y(Fg,ag,T-t,K)+(T-t)*ag).^2).*N6(d2(Fg,ag,T-t,K));
source_AA = (1/4)*rho*K*(T-t)^(3/2)*(ag.^2).*N3(d2(Fg,ag,T-t,K))...
    -(1/4)*rho*K*ag.*((y(Fg,ag,T-t,K)+(T-t)*ag).^2).*N4(d2(Fg,ag,T-t,K));
source3 = -(1/4)*rho*K*(T-t)^(1/2)*ag.*((-2/3)*(T-t)*ag+sqrt(T-t).*d2(Fg,ag,T-t,K)).*N3(d2(Fg,ag,T-t,K))...
    -(5/12)*rho*K*ag.*((y(Fg,ag,T-t,K)+(T-t)*ag).^2).*N4(d2(Fg,ag,T-t,K))...
    -(1/8)*(rho^3) * K *(T-t)^(1/2)*ag.*y(Fg,ag,T-t,K).*N5(d2(Fg,ag,T-t,K))...
    -(1/8)*(rho^3)*K*ag.*((y(Fg,ag,T-t,K)+(T-t)*ag).^2).*N6(d2(Fg,ag,T-t,K));
% source_A = (1/12)*K*(T-t)^(3/2)*N1(d2(Fg,ag,T-t,K)) ...
%     - (1/12)*K*(T-t)*(3*(T-t)*ag+sqrt(T-t).*d2(Fg,ag,T-t,K)).*N2(d2(Fg,ag,T-t,K))...
%     + K*(T-t)^(1/2)*((1/6)*((T-t)*ag+sqrt(T-t).*d2(Fg,ag,T-t,K)).^2+(1/8)*rho*rho*(T-t)).*N3(d2(Fg,ag,T-t,K))...
%     -(1/8)*(rho^2)*K*(T-t)*(y(Fg,ag,T-t,K)+2*(T-t)*ag).*N4(d2(Fg,ag,T-t,K))...
%     +(1/8)*(rho^2) * K *(T-t)^(1/2)*((T-t)*ag+sqrt(T-t).*d2(Fg,ag,T-t,K)).^2.*N5(d2(Fg,ag,T-t,K));

% %Check solutions
% BS = V_t + 0.5*a*a*(Fg.*Fg.*V_FF) + source;
% figure
% surf(Fg,tg,BS)
% title('BS1')
% 
% BS2 = W2_t + 0.5*a*a*(Fg.*Fg.*W2_FF) - source2;
% figure
% surf(Fg,tg,BS2)
% title('BS2')
% 
% %Check sources

error1 = source1 + rho*(ag.^2).*Fg.*W0_FA;
figure
surf(Fg,ag,error1)
title('error in source 1')

error2 = source2 + 0.5*(ag.^2).*W0_AA + rho*(ag.^2).*Fg.*W1_FA;
figure
surf(Fg,ag,error2)
title('error in source 2')

error_FA = source_FA + rho*(ag.^2).*Fg.*W2_FA  ;% + 0.5*(ag.^2).*W1_AA ;
figure
surf(Fg,ag,error_FA)
title('error in source ')

error_AA = source_AA + 0.5*(ag.^2).*W1_AA ;
figure
surf(Fg,ag,error_AA)
title('error in source ')

error3 = source3 + 0.5*(ag.^2).*W1_AA + rho*(ag.^2).*Fg.*W2_FA;
figure
surf(Fg,ag,error3)
title('error in source 3')


% BS = C_T-(0.5)*alpha*alpha*(S.^4).*C_SS;
% BS_K = C_T-(0.5)*alpha*alpha*((K).^2).*(K.^2).*C_KK;
% 
% Dupire = sqrt(2*C_T./((K.*K).*C_KK));
% Vidal = sqrt(2*C_T./((S.*S).*C_SS));
% LocalVol = (alpha)*K;
% figure
% plot(K,Dupire)
% hold on
% plot(K,LocalVol,'red')
% hold on
% plot(K,Vidal,'green')
% title(['S_0 = ' num2str(S0)])
% xlabel('K')
% ylabel('LocalVol')
% legend('Dupire formula','LocalVol','BlackScholes')
% figure
% plot(K,BS)
% 
% [MonteCarlo,confidence]  = C_MC_Dupire(100,1, 80 , alpha ,100000,1000)
% 
% Analytical = C(100,1,80)
% 
% [MonteCarlo,confidence]  = C_MC_Dupire(100,1, 100 , alpha ,100000,1000)
% 
% Analytical = C(100,1,100)
% 
% [MonteCarlo,confidence]  = C_MC_Dupire(100,1, 120 , alpha ,100000,1000)
% 
% Analytical = C(100,1,120)

