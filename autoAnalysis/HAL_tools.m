%% Collection of functions defining the HAL module
function handle = HAL_tools(functionName)
switch functionName
    case 'PCA'
        handle = @PCA;
    case 'eval'
        handle = @eval;
    case 'analyse'
        handle = @analyse;
    case 'analyseSampling'
        handle = @analyseSampling;
    case 'cubeSampling'
        handle = @cubeSampling;
    case 'filterOscil'
        handle = @filterOscil;
    case 'periodVsBeta'
        handle = @periodVsBeta;
    case 'equ'
        handle = @equ;
    case 'jac'
        handle = @jac;
    case 'testIntegrator'
        handle = @testIntegrator;
    case 'testSampling'
        handle = @testSampling;
    otherwise
        warning('A function named %s does not seem to exist',functionName)
        handle = 0;
end
end

%% Extract k's
function k = extractK(models)
k = zeros(length(models),length(models(1).k));
for ii = 1:length(models)
    k(ii,:) = models(ii).k;
end
end

%% PCA analysis
function [param,exponents,spearmanCorr] = PCA(modelSelection)
k = extractK(modelSelection);
per = arrayfun(@(model) model.period, modelSelection);

if sum(per<=0) > 0
    warning(['All data should have a period bigger then 0,'...
        ' did you remove non-oscillatory data points?'])
end

data = [k' ; per];

N = length(per);

% Moving analysis to the log-scale
data = log(data);

% Calculating the covariance matrix
covarianceMat = cov(data');

% Getting the smallest eigenvector
[V,D] = eig(covarianceMat);
exponents = D(1,1).*V(1:end-1,1);

% Calculating the predicting parameter (previously known as beta)
param = prod(k .^ (ones(N,1) * exponents'),2);

% Evaluating the predicting strenght of the resulting parameter
spearmanCorr = corr(param,per','type','Spearman');
end

%% Plot period vs beta
function periodVsBeta(models)
figure(1)
clf

[betaAnalytic,TAanalytic,TBanalytic] = analyticHAL(0.5:0.1:20);
plot(betaAnalytic,TAanalytic,'b-','linewidth',3)
hold on
plot(betaAnalytic,TBanalytic,'g-','linewidth',3)
plot(betaAnalytic,TAanalytic+TBanalytic,'k-','linewidth',3)

modelsSelection = filterOscil(models);
betas = arrayfun(@(model) model.beta, modelsSelection);
dm = arrayfun(@(model) model.k(2), modelsSelection);
phaseA = arrayfun(@(model) model.phaseA, modelsSelection);
phaseB = arrayfun(@(model) model.phaseB, modelsSelection);
periods = arrayfun(@(model) model.period, modelsSelection);
plot(betas,phaseA.*dm,'b.')
plot(betas,phaseB.*dm,'g.')
plot(betas,periods.*dm,'k.')

axis([0 10 0 20])

legend('length A phase','length B phase','period','length A phase','length B phase','period')

title(sprintf('(warnings=%i)',length(models)-length(modelsSelection)))
title('Matching of the different phases between numerical and analytical data')
xlabel('beta')
ylabel('period')

set(findall(gcf,'-property','FontSize'),'FontSize',25)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',10)
end

%% Analytical HAL
function [beta,TA,TB] = analyticHAL(TB)

equ = @(T1,T2) T1.*(1-exp(-T1-T2))./((1-exp(-T1)).*(-1+T2+exp(-T2))) + ...
    T2.*(1-exp(-T1-T2))./((T2-T2.^2./2).*(1-exp(-T1-T2))+(-exp(-T1)+T2.*exp(-T1)+1).*(exp(-T2)-1));

TA = zeros(size(TB));

for ii = 1:size(TB,2)
    T2local = TB(ii);
    localEqu = @(T1) equ(T1,T2local);
    
    TA(ii) = fzero(localEqu,10);
end

beta = TA.*(1-exp(-TA-TB))./((1-exp(-TA)).*(-1+TB+exp(-TB)));
end

%% Filter oscillators
function modelsSelection = filterOscil(models)
modelsSelection = models(arrayfun(@(model) model.oscil==1, models));
end

%% Analysis and sampling combined
function models = analyseSampling(baseModel,N,scale,simLength,filename)
models = cubeSampling(baseModel,N,scale);

parfor ii = 1:length(models)
    % assignement needed inside parfor eventhough model is handle (pointer)
    models(ii) = analyse(models(ii),simLength);
end
%arrayfun(@(model) analyse(model,simLength), models);

[path, ~, ~] = fileparts(mfilename('fullpath'));
save(strcat(path,'/data/',filename))
end

%% Analyse period of model
function model = analyse(model,simLength)

track = eval(model,simLength);
[model.numPer, model.phaseA, model.phaseB] = calcPeriod(track);

while model.numPer < 10 && track.x(end) < simLength*50
    track = odextend(track,[],track.x(end)+simLength);
    [model.numPer, model.phaseA, model.phaseB] = calcPeriod(track);
end

if model.numPer >= 10
    model.oscil = 1;
    % Only interesting if a limit cycle was found
    [model.maxA, model.minA, model.maxB, model.minB] = calcAmp(track);
else
    model.oscil = 0;
end
end

%% Calculate period
function [numPer, phaseA, phaseB] = calcPeriod(track)
% 0 => A phase; 1 => B phase
phase = track.y(3,:) < track.y(4,:);
% list of transition times from A phase to B phase
AtoB = track.x(diff(phase)==1);
% list of transition times from B phase to A phase
BtoA = track.x(diff(phase)==-1);

% if unequal periods of A and B remove one
if length(BtoA) < length(AtoB)
    AtoB(end) = [];
end

numPer = length(AtoB);

% break the evaluation if to little periods found
if(numPer<10)
    phaseA = 0; phaseB = 0;
    return
end

% estimate the length of the B phase by the last 5 instances
Bphases = BtoA-AtoB;
phaseB = mean(Bphases(end-5:end));

% estimate the length of the A phase by the last 5 instances
Aphases = AtoB(2:end)-BtoA(1:end-1);
phaseA = mean(Aphases(end-5:end));
end

%% Calculate amplitude of the limit cycle
function [maxA, minA, maxB, minB] = calcAmp(track)
% Calculate maxima and minima of the limit cycle for A and B
% (last 4/5-th to skip the transient)
maxA = max(track.y(3,floor(end*4/5):end));
minA = min(track.y(3,floor(end*4/5):end));
maxB = max(track.y(4,floor(end*4/5):end));
minB = min(track.y(4,floor(end*4/5):end));
end

%% Test sampling
function testSampling(baseModel,N,scale)

models = cubeSampling(baseModel,N,scale);

figure(999)
clf

for ii = 1:12
    p1 = arrayfun(@(model) model.k(ii          ), models);
    p2 = arrayfun(@(model) model.k(mod(ii,12)+1), models);
    
    subplot(3,4,ii)
    loglog(p1,p2,'b*')
    hold on
    loglog(baseModel.k(ii),baseModel.k(mod(ii,12)+1),'k*','markersize',10)
    
    xlabel(sprintf('param%i',ii))
    ylabel(sprintf('param%i',mod(ii,12)+1))
    axis tight
    
    % enforce 3 ticks
    set(gca,'XTick',logspace(log10(min(p1)*1.01),log10(max(p1)*0.99),3))
    set(gca,'YTick',logspace(log10(min(p2)*1.01),log10(max(p2)*0.99),3))
end

set(findall(gcf,'-property','FontSize'),'FontSize',15)
end

%% Algorithm to produce an Array of models in the neighbourhood of baseModel
function models = cubeSampling(baseModel,N,scale)
% Dimention of the space to be sampled
dim = length(baseModel.k);

% Random values between -1 and 1
r = rand(dim,N)*2-1;

% Params in a logaritmic cube aroud kInitial
params = exp(r*log(scale)+log(baseModel.k)'*ones(1,N))';

% Make models (without for-loop?)
models(N) = model();
for ii = N:-1:1
    models(ii).k = params(ii,:);
end
end

%% Test integrator
function testIntegrator(baseModel,simLength)
track = eval(baseModel,simLength);
plotTimeTracks(track);
end

%% Plot the time track of the concentrations
function plotTimeTracks(track)
figure(999)
clf
semilogy(track.x,track.y,'linewidth',5)
legend(equNames())
title('Concentration vs time')
xlabel('time')
ylabel('concentration')
set(findall(gcf,'-property','FontSize'),'FontSize',25)
end

%% Return the names of the different equations
function names = equNames()
names = {'G','M','A','B','AB'};
end

%% Integrate the model over time
function track = eval(baseModel,simLength)
% Use a stiff ode solver
track = ode23tb(...
    ... to integrate the HAL equation
    @(~,x)equ(x,baseModel.k), ...
    ... till tEnd and from x0
    [0 simLength], baseModel.x0, ...
    ... with access to the jacobian for performance
    odeset('Jacobian',@(~,x)jac(x,baseModel.k)));
end

%% ODE of the HAL module
function f = equ(x,k)
f = [k(1)*(1-x(1)) - k(7)*x(1)*x(3) ; ...
    k(3)*x(1) + k(8) * (1-x(1)) - k(2)*x(2) ; ...
    k(4)*x(2) - k(9)*x(3) - k(6)*x(3)*x(4) + k(12)*x(5) + k(1)*(1-x(1)) - k(7)*x(1)*x(3) ; ...
    k(5) - k(10) * x(4) - k(6)*x(3)*x(4) + k(12)*x(5) ; ...
    k(6)*x(3)*x(4) - k(12)*x(5) - k(11)*x(5)];
end

%% Jacobian of the HAL module
function J = jac(x,k)
J = [-k(1)-k(7)*x(3), 0    , -k(7)*x(1)               , 0               , 0           ; ...
    k(3)-k(8)       , -k(2), 0                        , 0               , 0           ; ...
    -k(1)-k(7)*x(3) , k(4) , -k(9)-k(6)*x(4)-k(7)*x(1), -k(6)*x(3)      , k(12)       ; ...
    0               , 0    , -k(6)*x(4)               , -k(10)-k(6)*x(3), k(12)       ; ...
    0               , 0    , k(6)*x(4)                , k(6)*x(3)       , -k(12)-k(11)];
end