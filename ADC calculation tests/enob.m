clear
s=csvread('v1.csv');
vin1 = s(:,2);
plot(vin1)
f1 = abs(fft(vin1));
plot(20*log10(f1))
psig1 = f1(2)^2 + f1(128)^2;
pt1 = sum(f1.^2);
snr1 = psig1/(pt1-psig1);
enob1 = (10*log10(snr1) - 1.76)/6.02;

%% T_test
s=csvread('v2.csv');
vin2 = s(:,2);
% plot(vin2)
f2 = abs(fft(vin2));
%% es
% plot(20*log10(f1))
psig2 = f2(64)^2 + f2(66)^2;
pt2 = sum(f2.^2);
snr2 = psig2/(pt2-psig2);
enob2 = (10*log10(snr2) - 1.76)/6.02;
