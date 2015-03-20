clear all
close all
clc

dxdt = @(~,x,k) ...
    [k(1)*(1-x(1)) - k(7)*x(1)*x(3) ; ...
    k(3)*x(1) + k(8) * (1-x(1)) - k(2)*x(2) ; ...
    k(4)*x(2) - k(9)*x(3) - k(6)*x(3)*x(4) + k(12)*x(5) + k(1)*(1-x(1)) - k(7)*x(1)*x(3) ; ...
    k(5) - k(10) * x(4) - k(6)*x(3)*x(4) + k(12)*x(5) ; ...
    k(6)*x(3)*x(4) - k(12)*x(5) - k(11)*x(5)];

params = [0.01, 0.05,  2,     25,  100,  1000, 0.01, 0.001, 0.001, 0.001, 1,   0.001];

[time, conc] = ode23tb(@(t,x)dxdt(t,x,params),[0 1000],[1 0 0 0 0]);

% remove transiant and scale time
time = time(round(end*3/4):end);

% shift and scale concentrations
conc = conc(round(end*3/4):end,1:5);
maxProt = max(max(conc(:,3:5)));
conc(:,3:5) = conc(:,3:5)./maxProt*0.9+1;
maxMrna = max(max(conc(:,2)));
conc(:,2) = conc(:,2)./maxMrna*0.4+0.5;
maxGene = max(max(conc(:,1)));
conc(:,1) = conc(:,1)./maxGene*0.4;

% plot the time tracks
figure(7)
clf

rectangle('Position', [time(1) 0.45 time(end)-time(1) 0.05],'FaceColor',[0.6 0.6 0.6])
hold on
rectangle('Position', [time(1) 0.95 time(end)-time(1) 0.05],'FaceColor',[0.6 0.6 0.6])
plot(time,conc,'linewidth',5)

% plot subphases
[~,phase_switch] = findpeaks(conc(:,3));
scaled_period = (phase_switch(end)-phase_switch(end-1))*params(4)
[~,phase_switch_extra] = findpeaks(conc(:,4));
phase_switch = [phase_switch; phase_switch_extra];
for ii = 1:length(phase_switch)
    plot([time(phase_switch) time(phase_switch)],[0 2],'k-')
end

% plot phases
[~,phase_switch] = findpeaks(min(conc(:,4),conc(:,3)));
for ii = 1:length(phase_switch)
    plot([time(phase_switch) time(phase_switch)],[0 2],'k-','linewidth',3)
end

set(gca,'ytick',[0 0.2 0.4 0.5 0.7 0.9 1 1.225 1.45 1.675 1.9],...
    'yticklabel',round([0 maxGene/2 maxGene 0 maxMrna/2 maxMrna 0 maxProt/4 maxProt/2 maxProt*3/4 maxProt]*100)/100,...
    'fontsize',20)
legend({'gene','mRNA','protein A','molecule B','complex AB'},'fontsize',30)

xlabel('time (min)','fontsize',30)
ylabel('concentration (molecule/\mu m^3)','fontsize',30)
title('Time track of HAL example','fontsize',30)
box on

axis tight

beta = params(5).*params(3).*params(1)./params(7)./params(4)