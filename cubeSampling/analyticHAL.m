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