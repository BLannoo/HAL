% make a random walk in the parameter space of HAL

addpath(genpath('~/Documents/Research/HAL_Enrico/'))

clear all
close all
clc

%% Settings

N = 1000; % Number of data points (~2sec per data point)
scale = 1.1;
seed = 2;

kInitial = [0.1, 0.06, 0.007, 50, 4, 100, 0.1, 0.001, 0.001, 0.001, 1, 0.01];

%% Pre allocations

params = zeros(N,12);
Bperiod = zeros(1,N);
Aperiod = zeros(1,N);
beta = zeros(1,N);
numPeaks = zeros(1,N);
Amax = zeros(1,N);
Amin = zeros(1,N);
Bmax = zeros(1,N);
Bmin = zeros(1,N);

%% Select Params

rng(seed)
params = exp((rand(length(kInitial),N)*2-1)*log(scale)+log(kInitial)'*ones(1,N))';

%% Do calculation

tic
parfor i1 = 1:N
    [Aperiod(i1),Bperiod(i1),numPeaks(i1),Amax(i1),Amin(i1),Bmax(i1),Bmin(i1),beta(i1)] = evaluateHAL(params(i1,:));
end
toc

numWarning = sum(numPeaks<10)

%% sampling

figure(4)
clf

for i1 = 1:12
    subplot(3,4,i1)
    loglog(params((numPeaks>9),i1),params((numPeaks>9),mod(i1,12)+1),'g.')
    hold on
    loglog(params((numPeaks<=9),i1),params((numPeaks<=9),mod(i1,12)+1),'r.')
    loglog(kInitial(i1),kInitial(mod(i1,12)+1),'b*','markersize',10)
    xlabel(sprintf('param%i',i1))
    ylabel(sprintf('param%i',mod(i1,12)+1))
    axis tight
end

set(findall(gcf,'-property','FontSize'),'FontSize',25)

%% Calc analytic phase lengths

[betaAnalytic,TAanalytic,TBanalytic] = analyticHAL(0.5:0.1:5);

%% Plot phase lengths

figure(1)
clf
hold all

semilogx(betaAnalytic,TAanalytic,'b-','linewidth',5)
semilogx(betaAnalytic,TBanalytic,'g-','linewidth',5)
semilogx(betaAnalytic,TAanalytic+TBanalytic,'k-','linewidth',5)

semilogx(beta(numPeaks>9),Aperiod(numPeaks>9).*params(numPeaks>9,2)','b.')
semilogx(beta(numPeaks>9),Bperiod(numPeaks>9).*params(numPeaks>9,2)','g.')
semilogx(beta(numPeaks>9),Aperiod(numPeaks>9).*params(numPeaks>9,2)'+Bperiod(numPeaks>9).*params(numPeaks>9,2)','k.')

axis([min(beta(numPeaks>9)) max(beta(numPeaks>9)) 0 5])

legend('length A phase','length B phase','period', ...
    'length A phase','length B phase','period')

title('Comparison phase lengths')
xlabel('beta')
ylabel('period')

set(findall(gcf,'-property','FontSize'),'FontSize',25)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',25)

%% Plot Amplitude distributions

figure(2)
clf
subplot(1,3,1)
boxplot([Amax(numPeaks>9);Bmax(numPeaks>9)]','labels',{'A','B'})
title('max value')
subplot(1,3,2)
boxplot([Amin(numPeaks>9);Bmin(numPeaks>9)]','labels',{'A','B'})
title('min value')
subplot(1,3,3)
boxplot([Amax(numPeaks>9)-Amin(numPeaks>9);Bmax(numPeaks>9)-Bmin(numPeaks>9)]','labels',{'A','B'})
title('abs amplitude')

set(findall(gcf,'-property','FontSize'),'FontSize',25)

%% Plot coefficient of variation

CV = @(x) std(x(numPeaks>9))./mean(x(numPeaks>9));

CV_Aperiod = CV(Aperiod);
CV_Bperiod = CV(Bperiod);
CV_beta = CV(beta);

CV_par = [];
for i1 = 1:12
    CV_par = [CV_par CV(params(:,i1))];
end

figure(3)
clf
boxplot([CV_Aperiod+zeros(size(CV_par));CV_Bperiod+zeros(size(CV_par));...
    CV_par;CV_beta+zeros(size(CV_par))]', ...
    'labels',{'length A phase','length B phase','params','beta'})

ylabel('CV')

set(findall(gcf,'-property','FontSize'),'FontSize',25)
