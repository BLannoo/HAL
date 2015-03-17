function [Aperiod,Bperiod,numPeaks,Amax,Amin,Bmax,Bmin,beta] = evaluateHAL(k)
    beta = k(1)*k(3)*k(4)/(k(5)*k(2)^2);
    
    track = getTracks(k,1000);
    [Aperiod,Bperiod,numPeaks] = getPeriods(track);
    
    if numPeaks < 10
        fprintf('Warning: <10 periods (%i) simlength*2\n',numPeaks)
        track = getTracks(k,2000);
        [Aperiod,Bperiod,numPeaks] = getPeriods(track);
        if numPeaks < 10
            fprintf('Warning: Still <10 periods (%i) returning\n',numPeaks)
            Amax=0;Amin=0;Bmax=0;Bmin=0;
            return
        end
    end

    Amax = max(track.y(3,floor(end*4/5):end));
    Amin = min(track.y(3,floor(end*4/5):end));
    Bmax = max(track.y(4,floor(end*4/5):end));
    Bmin = min(track.y(4,floor(end*4/5):end));
end

function track = getTracks(k,simLength)
    dxdt = @(~,x) ...
        [k(1)*(1-x(1)) - k(7)*x(1)*x(3) ; ...
        k(3)*x(1) + k(8) * (1-x(1)) - k(2)*x(2) ; ...
        k(4)*x(2) - k(9)*x(3) - k(6)*x(3)*x(4) + k(12)*x(5) + k(1)*(1-x(1)) - k(7)*x(1)*x(3) ; ...
        k(5) - k(10) * x(4) - k(6)*x(3)*x(4) + k(12)*x(5) ; ...
        k(6)*x(3)*x(4) - k(12)*x(5) - k(11)*x(5)];

    % starting A at 1 garanties A phase comes first
    track = ode23tb(dxdt,[0 simLength],[1 0 1 0 0]);
end

function [Aperiod,Bperiod,numPeaks] = getPeriods(track)
    phase = track.y(3,:) < track.y(4,:);
    AtoB = track.x(diff(phase)==1);
    BtoA = track.x(diff(phase)==-1);

    % if unequal periods of A and B remove one
    if length(BtoA) < length(AtoB)
        AtoB(end) = [];
    end
    
    numPeaks = length(AtoB);
    if(numPeaks<10)
        Aperiod=0;Bperiod=0;
        return
    end

    Bperiods = BtoA-AtoB;
    Bperiod = mean(Bperiods(end-5:end));
    
    Aperiods = AtoB(2:end)-BtoA(1:end-1);
    Aperiod = mean(Aperiods(end-5:end));
end