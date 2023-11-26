s=csvread('data19M.csv');
vout=s(:,2)/1.8;

s=abs(fft(vout));

psig=s(503)^2+s(523)^2;
ptotal=sum(s.^2)-s(1)^2;
figure
plot(20*log10(s/s(503)));
figure
plot(vout);
sndr=10*log10(psig/(ptotal-psig));
ENOB = (sndr-1.76)/6.02;

s1=s;
s1(1)=0;
s1(503)=0;
s1(523)=0;
sfdr=20*log10(s(503)/max(s1));