clear all
close all
clc
B=8;
W=B-1;
LSB=1/(200*(2^W));
R=2*(1/LSB)+1;
OUT=zeros(1,R);
OUT1=zeros(1,R);
OUT2=zeros(1,R);
OUT3=zeros(1,R);
OUT4=zeros(1,R);
OUT5=zeros(1,R);
OUT6=zeros(1,R);
bi1=zeros(1,R);
bi2=zeros(1,R);
bi3=zeros(1,R);
bi4=zeros(1,R);
bi5=zeros(1,R);
bi6=zeros(1,R);
bi7=zeros(1,R);
IN=-1:LSB:1;  
REF=1;
KH=0;
AC=50;
DNL=0;
DN=0;
DNLP=0;
INL=0;
SIGMA=0.03;
%Y=0;
Y=random('NORMAL',0,SIGMA);%Y = random(name,A,B) returns random numbers Y from a two-parameter family of distributions. Parameter values for the distribution are given in A and B.
P=zeros(1,256);
d=zeros(1,256);
INP=zeros(1,256);
for S=1:1:R
[OUT1(1,S), bi1(1,S)]=TOWHALF (IN(1,S)  , REF, 3, Y);
[OUT2(1,S), bi2(1,S)]=TOWHALF (OUT1(1,S), REF, 2, Y);
[OUT3(1,S), bi3(1,S)]=TOWHALF (OUT2(1,S), REF, 1, Y);
[bi7(1,S)]=FLASH2BIT (OUT3(1,S), REF);
OUT(1,S)= bi1(1,S) + bi2(1,S) + bi3(1,S)+ bi7(1,S);
end
OUTP=OUT;
figure
plot(IN,OUT1)
figure
plot(IN,OUT2)
figure
plot(IN,OUT3)
figure
plot(IN,OUTP)
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
DN=hist(OUT,code)/200-1;
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