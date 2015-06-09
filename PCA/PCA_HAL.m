%% PCA HAL cube sampling
clear all
close all
clc

% Getting the basic data (rates and period)
load('N1000_scale2.0_dataSet20.mat','params','Aperiod','Bperiod')
period = (Aperiod+Bperiod).*params(:,2)';
clear Aperiod Bperiod

% Eliminating non-oscillatory data points
k = params(period>0,:);
per = period(period>0);

[param,exponents,spearmanCorr] = GRN_tools.PCA(k,per);

% Checking the result
figure(2)
clf
plot(param,per,'.')
xlabel('P')
ylabel('period*\delta_m')
axis([min(param) max(param) 0 max(per)])




































