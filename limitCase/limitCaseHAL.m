%% Finding a limit case for the HAL model parameters

%% Settings

N = 9;
for par = 1:7; % 8-12 do not do anything

%        th      dm   km ka    kb  kab     al    kma  da   db   dab    diss
kBase = [0.0001, 0.1, 2, 1000, 25, 100000, 1000, 0.0, 0.0, 0.0, 10000, 0.0];

ks = ones(N,1)*kBase;
ks(:,par) = kBase(par) .* (2).^((1:N)-(N+1)/2);

for i1 = 1:N

    k = ks(i1,:);
    
    % Integrating
    dxdt = @(~,x) ...
        [k(1)*(1-x(1)) - k(7)*x(1)*x(3) ; ...
        k(3)*x(1) + k(8) * (1-x(1)) - k(2)*x(2) ; ...
        k(4)*x(2) - k(9)*x(3) - k(6)*x(3)*x(4) + k(12)*x(5) + k(1)*(1-x(1)) - k(7)*x(1)*x(3) ; ...
        k(5) - k(10) * x(4) - k(6)*x(3)*x(4) + k(12)*x(5) ; ...
        k(6)*x(3)*x(4) - k(12)*x(5) - k(11)*x(5)];

    % starting A at 1 garanties A phase comes first
    track = ode23tb(dxdt,[0 1000],[1 0 1 0 0]);

    % phase evaluation
    beta(i1) = k(1)*k(3)*k(4)/(k(5)*k(2)^2);

    phase = track.y(3,:) < track.y(4,:);
    AtoB = track.x(diff(phase)==1);
    BtoA = track.x(diff(phase)==-1);

    % if unequal periods of A and B remove one
    if length(BtoA) < length(AtoB)
        AtoB(end) = [];
    end
    
    numPeaks = length(AtoB);
    if(numPeaks<10)
        Aperiod(i1)=0;Bperiod(i1)=0;
        continue
    end

    Bperiods = BtoA-AtoB;
    Bperiod(i1) = mean(Bperiods(end-5:end));

    Aperiods = AtoB(2:end)-BtoA(1:end-1);
    Aperiod(i1) = mean(Aperiods(end-5:end));
    
end

% Plotting

figure(par)
clf

[betaAnalytic,TAanalytic,TBanalytic] = analyticHAL(0.5:0.1:10);
hold all

semilogx(betaAnalytic,TAanalytic,'b-')
semilogx(betaAnalytic,TBanalytic,'g-')
semilogx(betaAnalytic,TAanalytic+TBanalytic,'k-')

semilogx(beta,Aperiod.*ks(:,2)','bO')
semilogx(beta,Bperiod.*ks(:,2)','gO')
semilogx(beta,(Aperiod+Bperiod).*ks(:,2)','kO')

semilogx(beta,Aperiod.*ks(:,2)','bO')
semilogx(beta,Bperiod.*ks(:,2)','gO')
semilogx(beta,(Aperiod+Bperiod).*ks(:,2)','kO')

set(findall(gcf,'-property','FontSize'),'FontSize',15)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',15)
set(findall(gcf,'-property','linewidth'),'linewidth',5)


end


































