%% This code verifies that the correction terms W_1, W_2, and W_3 solve
% the corresponding PDE's. See stochasticVol.pdf for more info

%clear
clear all
close all
clc

%Parameters
a = 0.20 ;
T = 1 ;
rho = -0.10 ;
K = 30 ;
n = 2 ;


%functions
d2 = @(F,alpha,tau,strike) ( log( (1.0/strike)*F ) - 0.5*(alpha^2).*tau )./(alpha*sqrt(tau)) ;
d1 = @(F,alpha,tau,strike) ( log( (1.0/strike)*F ) + 0.5*(alpha^2).*tau )./(alpha*sqrt(tau)) ;
y = @(F,alpha,tau,strike) ( log( (1.0/strike)*F ) - 0.5*(alpha^2).*tau )./(alpha);

N = @(x) normcdf(x) ;
N1 = @(x) normpdf(x) ;
N2 = @(x) -x.*normpdf(x) ;
N3 = @(x) (x.^2 - 1.0 ).*normpdf(x) ;
N4 = @(x) (-x.^3+3*x).*normpdf(x) ;
N5 = @(x) (x.^4 - 6*x.^2 + 3).*normpdf(x) ;
N6 = @(x) (-x.^5 + 10*x.^3 - 15*x).*normpdf(x) ;



W0 = @(F,t) F.*N(d1(F,a,T-t,K)) - K*N(d2(F,a,T-t,K)) 
W1 = @(F,t) 0.5*rho*K*a*((T-t).*N2(d2(F,a,T-t,K)))
W2 = @(F,t) -(1/2)*a*K*((1/3)*a*((T-t).^2).*N2(d2(F,a,T-t,K))...
    + (1/3)*((T-t).^(3/2)).*d2(F,a,T-t,K).*N2(d2(F,a,T-t,K))...
    -(1/6)*((T-t).^(3/2)).*N1(d2(F,a,T-t,K)))...
    -(1/2)*rho*rho*a*K*((1/4)*a*((T-t).^2).*N4(d2(F,a,T-t,K))...
    + (1/4)*((T-t).^(3/2)).*d2(F,a,T-t,K).*N4(d2(F,a,T-t,K))...
    -(1/4)*((T-t).^(3/2)).*N3(d2(F,a,T-t,K)))
W3 = @(F,t) (1/24)*rho*K*a*((T-t).^2).*N2(d2(F,a,T-t,K))...
    -(1/16)*rho*K*a*(T-t).^(3/2).*(y(F,a,T-t,K)+(8/3)*a*(T-t)).*N3(d2(F,a,T-t,K))...
    +(1/12)*rho*K*a*(T-t).*((y(F,a,T-t,K)+a*(T-t)).^2+(1/5)*rho^2*(T-t)+(1/4)*(T-t)).*N4(d2(F,a,T-t,K))...
    -(1/60)*rho^3*K*a*(T-t).^(3/2).*(y(F,a,T-t,K)+(5/2)*a*(T-t)).*N5(d2(F,a,T-t,K))...
    +(1/48)*rho^3*K*a*(T-t).*((y(F,a,T-t,K)+a*(T-t)).^2+(1/5)*(T-t)).*N6(d2(F,a,T-t,K));
V = @(F,t) -((1/((n+2)*(n+3)))*(T-t) + (1/(n+3))*(y(F,a,T-t,K) + a*(T-t)).^2).*((T-t).^(n+1-2)).*N4(d2(F,a,T-t,K))...
    +2*3*((1/((n+2)*(n+3)))*(y(F,a,T-t,K) + a*(T-t))).*((T-t).^(n+1-3/2)).*N3(d2(F,a,T-t,K))...
    -2*3*2*(1/((n+1)*(n+2)*(n+3))).*((T-t).^(n)).*N2(d2(F,a,T-t,K))
V1 = @(F,t) -(1/(n+1))*((T-t).^(n-1)).*N4(d2(F,a,T-t,K))

%Plot
F = 20:1:60 ;
t = 0:0.05:.95 ;
[Fg, tg] = meshgrid(F,t) ;
z = W1(Fg,tg) ;
surf(Fg,tg,z)
title('W^1')

% Check source term
h = 0.00001;
dt = 0.000000001;
%W3_FF = (1/(h*h))*(W3(Fg+h,tg)-2*W3(Fg,tg)+W3(Fg-h,tg)) ;
%W3_t = (1.0/dt)*(W3(Fg,tg+dt)-W3(Fg,tg));
V_FF = (1/(h*h))*(V(Fg+h,tg)-2*V(Fg,tg)+V(Fg-h,tg)) ;
V_t = (1.0/dt)*(V(Fg,tg+dt)-V(Fg,tg));
W3_FF = (1/(h*h))*(W3(Fg+h,tg)-2*W3(Fg,tg)+W3(Fg-h,tg)) ;
W3_t = (1.0/dt)*(W3(Fg,tg+dt)-W3(Fg,tg));
W2_FF = (1/(h*h))*(W2(Fg+h,tg)-2*W2(Fg,tg)+W2(Fg-h,tg)) ;
W2_t = (1.0/dt)*(W2(Fg,tg+dt)-W2(Fg,tg));
W1_FF = (1/(h*h))*(W1(Fg+h,tg)-2*W1(Fg,tg)+W1(Fg-h,tg)) ;
W1_t = (1.0/dt)*(W1(Fg,tg+dt)-W1(Fg,tg));
W0_FF = (1/(h*h))*(W0(Fg+h,tg)-2*W0(Fg,tg)+W0(Fg-h,tg)) ;
W0_t = (1.0/dt)*(W0(Fg,tg+dt)-W0(Fg,tg));
sourceV1 = ((T-tg).^(n-2)).*N4(d2(Fg,a,T-tg,K));
source = ((T-tg).^(n-2)).*((a*(T-tg)+ y(Fg,a,T-tg,K)).^2).*N4(d2(Fg,a,T-tg,K));
source1 = -rho*a*K*N2(d2(Fg, a, T-tg, K)); 
source2 = (1/2)*a*K*(a*(T-tg)+sqrt(T-tg).*d2(Fg,a,T-tg,K)).*N2(d2(Fg,a,T-tg,K))...
    +(1/2)*rho*rho*a*K*(a*(T-tg)+sqrt(T-tg).*d2(Fg,a,T-tg,K)).*N4(d2(Fg,a,T-tg,K));
source3 = -(1/4)*rho*a*K*(T-tg).^(1/2).*((-2/3)*a*(T-tg)+sqrt(T-tg).*d2(Fg,a,T-tg,K)).*N3(d2(Fg,a,T-tg,K))...
    -(5/12)*rho*a*K*((y(Fg,a,T-tg,K)+(T-tg)*a).^2).*N4(d2(Fg,a,T-tg,K))...
    -(1/8)*(rho^3) * a * K *(T-tg).^(1/2).*y(Fg,a,T-tg,K).*N5(d2(Fg,a,T-tg,K))...
    -(1/8)*(rho^3)*K*a*((y(Fg,a,T-tg,K)+a*(T-tg)).^2).*N6(d2(Fg,a,T-tg,K));
%Check
BS = V_t + 0.5*a*a*(Fg.*Fg.*V_FF) - source;
figure
surf(Fg,tg,BS)
title('BS1')

BS2 = W2_t + 0.5*a*a*(Fg.*Fg.*W2_FF) - source2;
figure
surf(Fg,tg,BS2)
title('BS2')

BS3 = W3_t + 0.5*a*a*(Fg.*Fg.*W3_FF) - source3;
figure
surf(Fg,tg,BS3)
title('BS3')


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

