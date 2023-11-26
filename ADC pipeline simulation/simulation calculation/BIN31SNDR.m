s=csvread('data1.2M.csv');
vout=s(:,2)/1.8;

s=abs(fft(vout));

psig=s(32)^2+s(994)^2;
ptotal=sum(s.^2)-s(1)^2;
figure
plot(20*log10(s/s(32)));
figure
plot(vout);
sndr=10*log10(psig/(ptotal-psig));
ENOB = (sndr-1.76)/6.02;

s1=s;
s1(1)=0;
s1(32)=0;
s1(994)=0;
sfdr=20*log10(s(32)/max(s1));
        
