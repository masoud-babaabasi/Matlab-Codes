s=csvread('ramp.csv');
vout=s(:,2);
vout=vout';
%figure
%plot(vout);

B=8;
W=B-1;
LSB=1/(200*(2^W));
R=length(vout);
 
REF=1;
KH=0;
AC=50;
DNL=0;
DN=0;
DNLP=0;
INL=0;

P=zeros(1,256);
d=zeros(1,256);
INP=zeros(1,256);

P=zeros(1,256);
OUTP=vout;
figure
plot(OUTP)

 for F=1:1:R
     for H=0:1:255
         if OUTP(1,R)==H
             P(1,H+1)=P(1,H+1)+1;
             INP(1,H+1)=IN(1,F);
         end
     end
 end
 for L=1:1:256
     if P(1,L)<0.5*AC  || P(1,L)>1.5*AC
         DNL=1;
     end
     if P(1,L)==0.5*AC || P(1,L)==1.5*AC
         KH=1;
     end
     if INP(1,L)>-1+(L+1)/128
         INL=1;
     end
 end
 DNLP=DNL;
 ERROR=KH;
 INLP=INL;
code=0:255;
DN=hist(vout,code)/200-1;
figure
plot(code,DN)
xlabel('Code'), ylabel('DNL');
title('DNL');
for i=3:1:256
    d(1,i)=0;
    for j=2:1:i-1
        d(1,i)= DN(1,j)+d(1,i);
    end
end
figure  
plot(d) 
xlabel('Code'), ylabel('INL');
title('INL');