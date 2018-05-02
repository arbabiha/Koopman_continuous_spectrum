function Spectrum_SP()
% data-driven approximation of Koopman continuous spectrum
% for "Study of dynamics in post-transient flows using Koopman mode
% decomposition", H. Arbabi & I. Mezic, Physical Review Fluids, 2017
% Hassan Arbabi, April 2018 arbabiha@gmail.com

addpath('./tools')

%% data generation
disp('generating data ....')
% Arnold's cat map
X=CatMapData(1e5);

% Lorenz system
dt = .2;   % sampling interval
fs = 1/dt;  % sampling frequency
Y=LorenzData(dt,5e4);


% cavity flow at Re=30k
load('cavity_obs.mat')
Z = obs_raw_Re30k;
ws= 2*pi*10;    % sampling frequency
%% analysis
disp('computing the spectrum ...')
% Welsch Method 
M = 80;     % length of each subsample
v = blackman(M);    % the window function
K = M/2;    % overlap length
L = M*4;   % fft grid

phi1 = welch_estimator(X(:,1),v,K,L);
phi2 = bartlett_estimator(X(:,2),M,L);
phi3 = welch_estimator(Y(:,1),v,K,L);
phi4 = welch_estimator(Z(:),v,K,L);


%% plots & comparison
set(0,'defaultTextInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex');


% analytical solution for cat map densities
theta = linspace(0,2*pi,100);
r1 = 5/4 + cos(theta);
r2 = 21/16 + (5/4)*cos(theta) + (1/2)*cos(2*theta);


figure(),clf
subplot(2,2,1)
plot(theta,r1,'k','LineWidth',2),hold on
plot(linspace(0,2*pi,L),phi1,'--','LineWidth',2)
xlim([0,2*pi]), xlabel('$\theta$'), legend('analytical','Welch')
title('cat map, $\rho(g_1)$')
subplot(2,2,2)
plot(theta,r2,'k','LineWidth',2),hold on
plot(linspace(0,2*pi,L),phi2,'--','LineWidth',2)
xlim([0,2*pi]), xlabel('$\theta$'), legend('analytical','Bartlett')
title('cat map, $\rho(g_2)$')


subplot(2,2,3)
semilogy(linspace(0,2*pi*fs,L),phi3/fs,'LineWidth',2) % we have to scale according to sampling frequency
xlim([0 fs*pi])
title('Lorenz, $\rho(x_1)$'), xlabel('$\omega$')
legend('Welch')


subplot(2,2,4)
semilogy(linspace(0,ws,L),phi4*2*pi/ws,'LineWidth',2) % we have to scale according to sampling frequency
xlim([0 ws/2])
title('cavity random observable, $\rho(\psi_1)$'),xlabel('$\omega$')
legend('Welch')


end


function [Y]=LorenzData(dt,n)
%% creates samples of data from Lorenz chaotic attractor
% dt : sampling interval 
% n  : number of samples

% Lorenz chaotic model 1963
sigma=10;
rho = 28;
beta = 8/3;

Lorenz = @(x) [sigma*(x(2)-x(1)); ...
               x(1)*(rho-x(3))-x(2);...
               x(1)*x(2)-beta*x(3)];

tspan = 0:dt:(n*dt + 20);   %% allow 20 seconds of transients


x0 = [0.1;0;0.1];
[tspan,Y]= ode45(@(t,y)Lorenz(y),tspan,x0);

Y = Y(tspan>20,:);
end


function [Y]=CatMapData(n)
%% creates samples of data from Arnold's cat map
% n  : number of samples

X = zeros(2,n);

X(:,1)=[.2;.37];

for i=1:n-1
    X(:,i+1)=mod([2 1; 1 1]*X(:,i),1);
end


% two observables
f1 = @(X) exp(2*pi*1i*(2*X(1,:)+X(2,:))) + (.5)*exp(2*pi*1i*(5*X(1,:)+3*X(2,:))) ;
f2 = @(X) f1(X)+ 0.25*exp(2*pi*1i*(13*X(1,:)+8*X(2,:)));
Y=[f1(X)',f2(X)'];
end