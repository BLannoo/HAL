%% This is the main script for the automated analysis of HAL

clear all
close all

%% Settings

N = 10000; % Number of data points
scale = 2; % Size of the sampled neighbourhood
seed = 2;
simLength = 100; % Initial simulation length

name={'th'  'dm'   'km'   'ka' 'kb'  'kab' 'al'  'kma'  'da'   'db'  'dab' 'diss'};
k = [[0.1,  0.06,  0.007, 50,  4,     100, 0.1,  0.001, 0.001, 0.001, 1,   0.01];
     [1,    0.36,  100,   100, 25459, 100, 1,    2,     0.08,  0.28,  1,   0.01];
     [0.1,  0.11,  48,    8.5, 1890,  100, 0.01, 5,     2.1,   0.92,  100, 100];
     [0.05, 0.33,  37,    1.4, 40,    100, 0.32, 0.9,   2.7,   0.04,  2.8, 14];
     [0.14, 0.4,   100,   0.9, 48,    4.9, 0.05, 3.3,   0.08,  0.16,  100, 100];
     [0.06, 0.13,  7.4,   1.9, 54,    17,  0.01, 1.8,   0.07,  0.26,  100, 100];
     [0.01, 0.1,   7.2,   0.2, 4.6,   8.8, 0.1,  0.22,  0.34,  0.02,  100, 100]; ... seed 139999
     [0.08, 0.04,  23,    0.5, 200,   0.8, 0.01, 7.27,  0.03,  0.02,  100, 100]; ... seed 130001
     [0.01, 0.05,  2,     25,  100,  1000, 0.01, 0.001, 0.001, 0.001, 1,   0.001]; ... analytic limit (close to biology)
     ];


allExponents = zeros(size(k));
for dataSet = 1:size(k,1);

%% Tests

% f = HAL_tools('testIntegrator'); f(baseModel,300);
% f = HAL_tools('testSampling'); f(baseModel,1000,scale);

%% Analysis

filename = sprintf('N%i_scale%.1f_dataSet%iv2.mat',N,scale,dataSet);

% Initial model
baseModel = model();
baseModel.k = k(dataSet,:);

if exist(filename,'file')
    load(filename)
else
    tic
    rng(seed)
    analyseSampling = HAL_tools('analyseSampling');
    models = analyseSampling(baseModel,N,scale,simLength,filename);
    toc
end

%% Plotting

f = HAL_tools('periodVsBeta');
f(models);

%% PCA

filterOscil = HAL_tools('filterOscil');
modelSelection = filterOscil(models);

PCA = HAL_tools('PCA');
[param,allExponents(dataSet,:),spearmanCorr] = PCA(modelSelection);


per = arrayfun(@(model) model.period, modelSelection);

% Checking the result
figure(100+dataSet)
clf
plot(param,per,'.')
title(sprintf('SpearmanCorr=%f',spearmanCorr))
xlabel('P')
ylabel('period*\delta_m')
axis([min(param) max(param) 0 max(per)])
set(findall(gcf,'-property','FontSize'),'FontSize',25)

end

%%

select = [1:12];

tmp = allExponents ./ (sign(allExponents(:,1)) * ones(1,12));
n = tmp ./ (max(abs(tmp),[],2) * ones(1,12));
figure(999)
clf
hist(n(:,select))
legend(name{select})













