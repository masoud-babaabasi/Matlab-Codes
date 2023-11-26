clear
s=csvread('Sample1024.csv');
vin1 = s(1:1024,2);
plot(vin1)
f1 = abs(fft(vin1));
psig1 = f1(505)^2 + f1(521)^2;
pt1 = sum(f1.^2);
snr1 = psig1/(pt1-psig1);
enob1 = (10*log10(snr1) - 1.76)/6.02;

v2 = vin1(1:128);
f2 = abs(fft(v2));
psig2 = f2(64)^2 + f2(66)^2;
pt2 = sum(f2.^2);
snr2 = psig2/(pt2-psig2);
enob2 = (10*log10(snr2) - 1.76)/6.02;