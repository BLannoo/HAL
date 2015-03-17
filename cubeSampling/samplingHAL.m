function [params,Aperiod,Bperiod,numPeaks,Amax,Amin,Bmax,Bmin,beta,numWarning] = samplingHAL(N,kInitial,scale,seed,fileName)

%% Pre allocations

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

parfor i1 = 1:N
    [Aperiod(i1),Bperiod(i1),numPeaks(i1),Amax(i1),Amin(i1),Bmax(i1),Bmin(i1),beta(i1)] = evaluateHAL(params(i1,:));
end

numWarning = sum(numPeaks<10);

save(filename)
