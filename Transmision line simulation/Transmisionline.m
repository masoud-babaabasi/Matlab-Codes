% This code is written to simulate the behavior of a ideal transmision line
%author: Masoud Babaabasi April 2021
clear
Zs = 50;
Zl = 100000;
Z0 = 32;
Vs = 1;
T_delay = 1e-9;
%stady state error 
error = 0.10; % 10 percent error
%number of points
np = 200;
%sampling time
ts = T_delay / np;
%*****
Vdc = Zl/(Zl+Zs)*Vs;
V1= Z0 / (Z0+Zs)*Vs;
K1 = (Zs - Z0)/(Zs + Z0);
K2 = (Zl - Z0)/(Zl + Z0);
%********
t=0;
ti = 0:ts:ts;
Vline = zeros(1,np);
Vr = zeros(1,np);
Vf = zeros(1,np);
VL = zeros(1,2);
Vi = zeros(1,2);
L=linspace(0,1,np);

figure('Name','transmition line')
subplot(3,1,1)
h1 = plot(L,Vline);
xticks([0 1])
xticklabels({'x=0','x=L'})
ylim([0 1.2*Vs])
ylabel('volts')
% title(sprintf('t=%f ns',t*1e9));
tick = {'t=0'};
subplot(3,1,2)
h2 = plot(ti,VL);
xticks(0);
xticklabels(tick)
ylim([0 1.2*Vs])
ylabel('V load')
subplot(3,1,3)
h3 = plot(ti,Vi);
xticks(0);
xticklabels(tick)
ylim([0 1.2*Vs])
ylabel('V source')
con=0;
Vf(1) = V1;
while con==0 || t ==0
    Vf(2:end) = Vf(1:end-1);
    Vr(1:end-1) = Vr(2:end);
    Vf(1) = V1 + K1*Vr(1);
    Vr(end) = Vf(end)*K2;
    Vline = Vf + Vr;
    VL(round(t/ts)+1) = Vline(end);
    Vi(round(t/ts)+1) = Vline(1);
    ti(round(t/ts)+1) = t;
    h1.YData = Vline;
    h2.YData = VL;
    h2.XData = ti;
    h3.YData = Vi;
    h3.XData = ti;
    subplot(3,1,1)
    ylim([min(0,1.2*min(Vline)) max([1.2*max(Vline) 1.2*Vs])])
    subplot(3,1,2)
    ylim([min(0,1.2*min(VL)) max([1.2*max(VL) 1.2*Vs])])
    xlim([0 ti(end)])
    subplot(3,1,3)
    ylim([min(0,1.2*min(Vi)) max([1.2*max(Vi) 1.2*Vs])])
    xlim([0 ti(end)])
    if mod(length(ti),np) == 0
        i = length(ti)/np;
        tick{i+1} = sprintf('t=%dtd',i);
        subplot(3,1,2)
        xticks(0:T_delay:i*T_delay);
        xticklabels(tick)
        subplot(3,1,3)
        xticks(0:T_delay:i*T_delay);
        xticklabels(tick)
    end
%     title(sprintf('t=%f ns',t*1e9));
    if abs(var(Vline)) < 1e-10
        if abs((Vline(1) - Vdc)/Vdc) < error
            con = 1;
        end
    end
    t = t + ts;
    pause(0.01);
end

