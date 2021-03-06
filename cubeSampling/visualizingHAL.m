% make a random walk in the parameter space of HAL

clear all
close all

%% Settings

N = 100; % Number of data points (~2sec per data point)
scale = 2;
seed = 2;

%     th    dm     km     ka   kb     kab  al    kma    da     db     dab  diss
k = [[0.1,  0.06,  0.007, 50,  4,     100, 0.1,  0.001, 0.001, 0.001, 1,   0.01];
     [1,    0.36,  100,   100, 25459, 100, 1,    2,     0.08,  0.28,  1,   0.01];
     [0.1,  0.11,  48,    8.5, 1890,  100, 0.01, 5,     2.1,   0.92,  100, 100];
     [0.05, 0.33,  37,    1.4, 40,    100, 0.32, 0.9,   2.7,   0.04,  2.8, 14];
     [0.14, 0.4,   100,   0.9, 48,    4.9, 0.05, 3.3,   0.08,  0.16,  100, 100];
     [0.06, 0.13,  7.4,   1.9, 54,    17,  0.01, 1.8,   0.07,  0.26,  100, 100];
     [0.01, 0.1,   7.2,   0.2, 4.6,   8.8, 0.1,  0.22,  0.34,  0.02,  100, 100]; ... seed 139999
     [0.08, 0.04,  23,    0.5, 200,   0.8, 0.01, 7.27,  0.03,  0.02,  100, 100]; ... seed 130001
     [0.0001, 0.1, 2,    1000, 25, 100000, 1000, 0.0,   0.0,   0.0, 10000, 0.0]; ... analytic limit
     [0.05, 0.06,  2,     50,  100,   100, 0.1,  0,     0,     0,     1,   0]; ... analytic limit (close to biology)
     [0.05, 0.2,   2,     50,  100,   100, 0.05, 0,     0,     0,     1,   0]; ... analytic limit (close to biology)
     [0.01, 0.1,   2,     50,  25,    100, 0.05, 0,     0,     0,     1,   0]; ... analytic limit (close to biology)
     [0.01, 0.1,   2,     50,  50,   1000, 0.05, 0,     0,     0,     1,   0]; ... analytic limit (close to biology)
     [0.01, 0.05,  2,     50,  100,  1000, 0.05, 0,     0,     0,     1,   0]; ... analytic limit (close to biology)
     [0.01, 0.05,  2,     50,  100,  1000, 0.05, 0,     0,     0,     1,   0]; ... analytic limit (close to biology)
     [0.01, 0.05,  2,     25,  100,  1000, 0.05, 0,     0,     0,     1,   0]; ... analytic limit (close to biology)
     [0.01, 0.075, 2,     25,  100,  1000, 0.05, 0,     0,     0,     1,   0]; ... analytic limit (close to biology)
     [0.01, 0.075, 2,     25,  100,  1000, 0.01, 0,     0,     0,     1,   0]; ... analytic limit (close to biology)
     [0.01, 0.05,  2,     25,  100,  1000, 0.01, 0,     0,     0,     1,   0]; ... analytic limit (close to biology)
     [0.01, 0.05,  2,     25,  100,  1000, 0.01, 0.001, 0.001, 0.001, 1,   0.001]; ... analytic limit (close to biology)
     ];

%% Get data

dataSet = 20;

filename = sprintf('N%i_scale%.1f_dataSet%i.mat',N,scale,dataSet);

kInitial = k(dataSet,:);

% if exist(filename,'file')
%     load(filename)
% else
    tic
    [params,Aperiod,Bperiod,numPeaks,Amax,Amin,Bmax,Bmin,beta,numWarning] = samplingHAL(N,kInitial,scale,seed,filename);
    toc
% end

%% sampling

% figure(100+dataSet)
% clf
% 
% for i1 = 1:12
%     subplot(3,4,i1)
%     loglog(params((numPeaks>9),i1),params((numPeaks>9),mod(i1,12)+1),'g.')
%     hold on
%     loglog(params((numPeaks<=9),i1),params((numPeaks<=9),mod(i1,12)+1),'r.')
%     loglog(kInitial(i1),kInitial(mod(i1,12)+1),'b*','markersize',10)
%     xlabel(sprintf('param%i',i1))
%     ylabel(sprintf('param%i',mod(i1,12)+1))
%     axis tight
% end
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',25)

%% Plot phase lengths

[betaAnalytic,TAanalytic,TBanalytic] = analyticHAL(0.5:0.1:20);

figure(1)
clf

plot(betaAnalytic,TAanalytic,'b-','linewidth',3)
hold all
plot(betaAnalytic,TBanalytic,'g-','linewidth',3)
plot(betaAnalytic,TAanalytic+TBanalytic,'k-','linewidth',3)

plot(beta(numPeaks>9),Aperiod(numPeaks>9).*params(numPeaks>9,2)','b.')
plot(beta(numPeaks>9),Bperiod(numPeaks>9).*params(numPeaks>9,2)','g.')
plot(beta(numPeaks>9),(Aperiod(numPeaks>9)+Bperiod(numPeaks>9)).*params(numPeaks>9,2)','k.')

axis([0 10 0 20])

legend('length A phase','length B phase','period','length A phase','length B phase','period')

title(sprintf('(warnings=%i)',numWarning))
title('Matching of the different phases between numerical and analytical data')
xlabel('beta')
ylabel('period')

set(findall(gcf,'-property','FontSize'),'FontSize',25)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',10)

%% Plot Amplitude distributions

% figure(2)
% clf
% subplot(1,3,1)
% boxplot([Amax(numPeaks>9);Bmax(numPeaks>9)]','labels',{'A','B'})
% title('max value')
% subplot(1,3,2)
% boxplot([Amin(numPeaks>9);Bmin(numPeaks>9)]','labels',{'A','B'})
% title('min value')
% subplot(1,3,3)
% boxplot([Amax(numPeaks>9)-Amin(numPeaks>9);Bmax(numPeaks>9)-Bmin(numPeaks>9)]','labels',{'A','B'})
% title('abs amplitude')
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',25)

%% Plot coefficient of variation

% CV = @(x) std(x(numPeaks>9))./mean(x(numPeaks>9));
% 
% CV_Aperiod = CV(Aperiod);
% CV_Bperiod = CV(Bperiod);
% CV_period = CV(Aperiod+Bperiod);
% CV_beta = CV(beta);
% 
% CV_par = [];
% for i1 = 1:12
%     CV_par = [CV_par CV(params(:,i1))];
% end
% 
% figure(3)
% %subplot(3,3,i2)
% boxplot([CV_Aperiod+zeros(size(CV_par));CV_Bperiod+zeros(size(CV_par));...
%     CV_period+zeros(size(CV_par));CV_par;CV_beta+zeros(size(CV_par))]' ...
%     ... ,'labels',{'length A phase','length B phase','period','params','beta'} ...
%     )
% 
% ylabel('CV')
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',25)
