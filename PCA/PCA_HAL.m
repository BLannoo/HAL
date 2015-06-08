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
data = [k' ; per];
clear period params

N = length(per);

% moving analysis to the log-scale
data = log(data);

% calculating the covariance matrix
covarianceMat = cov(data');

% Getting the smallest eigenvector
[V,D] = eig(covarianceMat);
n = D(1,1).*V(1:end-1,1)

% Calculating the predicting parameter (previously known as beta)
p = prod(k .^ (ones(N,1) * n'),2);

% Checking the result
figure(1)
clf
plot(p,per,'.')
xlabel('P')
ylabel('period*\delta_m')
axis([min(p) max(p) 0 max(per)])