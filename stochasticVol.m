%% This code verifies that the correction terms W_1 and W_2 solve
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


%functions
d2 = @(F,alpha,tau,strike) ( log( (1.0/strike)*F ) - 0.5*(alpha^2).*tau )./(alpha*sqrt(tau)) ;
d1 = @(F,alpha,tau,strike) ( log( (1.0/strike)*F ) + 0.5*(alpha^2).*tau )./(alpha*sqrt(tau)) ;

N = @(x) normcdf(x) ;
N1 = @(x) normpdf(x) ;
N2 = @(x) -x.*normpdf(x) ;
N3 = @(x) (x.^2 - 1.0 ).*normpdf(x) ;
N4 = @(x) (-x.^3+3*x).*normpdf(x) ;

W0 = @(F,t) F.*N(d1(F,a,T-t,K)) - K*N(d2(F,a,T-t,K)) 
W1 = @(F,t) 0.5*rho*K*a*((T-t).*N2(d2(F,a,T-t,K)))
W2 = @(F,t) -(1/2)*a*K*(-(1/3)*a*((T-t).^2).*N2(d2(F,a,T-t,K))...
    + (1/3)*((T-t).^(3/2)).*d2(F,a,T-t,K).*N2(d2(F,a,T-t,K))...
    -(1/6)*((T-t).^(3/2)).*N1(d2(F,a,T-t,K)))...
    +(1/2)*rho*rho*a*K*(-(1/4)*a*((T-t).^2).*N4(d2(F,a,T-t,K))...
    + (1/4)*((T-t).^(3/2)).*d2(F,a,T-t,K).*N4(d2(F,a,T-t,K))...
    -(1/4)*((T-t).^(3/2)).*N3(d2(F,a,T-t,K)))
%Plot
F = 20:1:60 ;
t = 0:0.05:.95 ;
[Fg, tg] = meshgrid(F,t) ;
y = W1(Fg,tg) ;
surf(Fg,tg,y)
title('W^1')

% Check source term
h = 0.00001;
dt = 0.000000001;
%W3_FF = (1/(h*h))*(W3(Fg+h,tg)-2*W3(Fg,tg)+W3(Fg-h,tg)) ;
%W3_t = (1.0/dt)*(W3(Fg,tg+dt)-W3(Fg,tg));
W2_FF = (1/(h*h))*(W2(Fg+h,tg)-2*W2(Fg,tg)+W2(Fg-h,tg)) ;
W2_t = (1.0/dt)*(W2(Fg,tg+dt)-W2(Fg,tg));
W1_FF = (1/(h*h))*(W1(Fg+h,tg)-2*W1(Fg,tg)+W1(Fg-h,tg)) ;
W1_t = (1.0/dt)*(W1(Fg,tg+dt)-W1(Fg,tg));
W0_FF = (1/(h*h))*(W0(Fg+h,tg)-2*W0(Fg,tg)+W0(Fg-h,tg)) ;
W0_t = (1.0/dt)*(W0(Fg,tg+dt)-W0(Fg,tg));
source = -rho*a*K*N2(d2(Fg, a, T-tg, K)); 
source2 = -(1/2)*a*K*(a*(T-tg)-sqrt(T-tg).*d2(Fg,a,T-tg,K)).*N2(d2(Fg,a,T-tg,K))...
    +(1/2)*rho*rho*a*K*(a*(T-tg)-sqrt(T-tg).*d2(Fg,a,T-tg,K)).*N4(d2(Fg,a,T-tg,K))
%-0.5*a*K*(a*(T-tg)-d2(Fg,a,T-tg,K).*sqrt(T-tg)).*N2(d2(Fg,a,T-tg,K))
 % +0.5*rho*rho*a*K*(a*(T-tg)-d2(Fg,a,T-tg,K).*sqrt(T-tg)).*N4(d2(Fg,a,T-tg,K));

%Check
BS = W1_t + 0.5*a*a*(Fg.*Fg.*W1_FF) - source;
figure
surf(Fg,tg,BS)
title('BS')

BS2 = W2_t + 0.5*a*a*(Fg.*Fg.*W2_FF) - source2;
figure
surf(Fg,tg,BS2)
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

