%% PCA HAL cube sampling

clear all
close all
clc

load('N1000_scale2.0_dataSet20.mat','params','Aperiod','Bperiod')
period = (Aperiod+Bperiod).*params(:,2)';
clear Aperiod Bperiod

k = params(period>0,:);
beta = k(:,1).*k(:,3).*k(:,4)./(k(:,5).*k(:,2).^2);
per = period(period>0);
clear period params

N = length(per);

%% Permute the measured periods randomly to break the correlation
% To test the algorithm on a similar non-correlated data set

%per = per(randperm(N));

%%

data = log([k' ; per]);
data = data - mean(data,2)*ones(1,N);
data = data ./ (std(data,0,2)*ones(1,N));
[eigVec, ~, eigVal] = princomp(data');

n = eigVal(1).*eigVec(:,1);

p = prod(k .^ (ones(N,1) * n'),2);

figure(1)
clf
loglog(p,per,'.')
corr(p,per')

figure(2)
clf
loglog(beta.^-1,per,'.')
corr(beta.^-1,per')