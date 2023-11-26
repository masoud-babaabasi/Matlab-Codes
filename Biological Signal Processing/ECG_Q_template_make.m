%Preprocessing
clear 
load Data.mat
N = length(data.signal);
Fs = data.samplingfreq;
ecg = data.signal;
%-----------------------------
t = 0:1/Fs:(N-1)/Fs;
f = 0:Fs/N:(N-1)*Fs/N;
figure('Name','Input Signal','NumberTitle','off'); 
plot(t,ecg)
%input spectrum---------------
Spectrum = fft(ecg);
figure('Name','Spectrm','NumberTitle','off'); 
subplot(2,1,1);
plot(f(1:N/2),abs(Spectrum(1:N/2))/N*2);
subplot(2,1,2);
plot(f(1:N/2),180/pi*phase(Spectrum(1:N/2)));
%Notch filter------------------
theta = 2 * pi * 50 / Fs;
r = 0.995;
b = [1 -2*cos(theta) 1];
a = [1 -2*r*cos(theta) r^2];
ecg2 = filter(b , a ,ecg);
figure('Name','Filteredsignal(Notch filter 50Hz)','NumberTitle','off'); 
plot(t,ecg2)
Spectrum2 = fft(ecg2);
figure('Name','Spectrm2(Notch filter 50Hz)','NumberTitle','off'); 
subplot(2,1,1);
plot(f(1:N/2),abs(Spectrum2(1:N/2))/N*2);
subplot(2,1,2);
plot(f(1:N/2),180/pi*phase(Spectrum2(1:N/2)));
%Highpass filter----------------
fcH=0.24;
i=floor(fcH*N/Fs)+1;
Spectrum3 = Spectrum2;
Spectrum3(1:i)=0;
Spectrum3(N-i+1:N)=0;
ecg3 = real(ifft(Spectrum3));
figure('Name','Spectrm3(Highpass filter)','NumberTitle','off'); 
subplot(2,1,1);
plot(f(1:N/2),abs(Spectrum3(1:N/2))/N*2);
subplot(2,1,2);
plot(f(1:N/2),180/pi*phase(Spectrum3(1:N/2)));
figure('Name','HighPass','NumberTitle','off'); 
plot(t,ecg3)
%Lowpass filter-----------------
fcL=35;
i=floor(fcL*N/Fs)+1;
Spectrum4 = Spectrum3;
Spectrum4(i:N-i+1)=0;
ecg4 = real(ifft(Spectrum4));
figure('Name','Spectrm4(Lowpass filter)','NumberTitle','off'); 
subplot(2,1,1);
plot(f(1:N/2),abs(Spectrum4(1:N/2))/N*2);
subplot(2,1,2);
plot(f(1:N/2),180/pi*phase(Spectrum4(1:N/2)));
figure('Name','LowPass','NumberTitle','off'); 
plot(t,ecg4)

%***********************************************************************
%feature extraction
pwave = zeros(N,1);
N_pw=40;
index=zeros(N_pw,9);

Pstart = [385 840 1270 1715 2181  2638 3076 3543 4017];
for i=1:9
    index(:,i) = Pstart(i):Pstart(i)+N_pw-1;
    pwave(index(:,i)) = ecg4(index(:,i));
end
%P wave mnually extracted
figure('Name','P waves','NumberTitle','off'); 
plot(t,pwave);

Temp = zeros(N_pw,1);
for i=1:9
    Temp = Temp + pwave(index(:,i)) ;
end
Temp = Temp/9;
figure('Name','Template','NumberTitle','off'); 
plot(t(1:N_pw),Temp);

std = 0.1;
noise = std*randn(N,1);
ecg_noise = ecg4 + noise;
[r,lags] = xcorr(ecg_noise,Temp);
figure('Name','ECG noise','NumberTitle','off'); 
plot(t,ecg_noise)
figure('Name','Cross correlation','NumberTitle','off'); 
plot(lags/Fs,r);

V_threshold=0.15;
hold on
plot(lags/Fs,V_threshold*ones(2*N-1,1))

p_pos = zeros(N,1);

%Q wave extraction -------------------
Qwave = zeros(N,1);
N_Qw=10;
indexQ=zeros(N_Qw,9);

Qstart = [435 892 1317 1765 2232  2687 3128 3593 4067];
for i=1:9
    indexQ(:,i) = Qstart(i):Qstart(i)+N_Qw-1;
    Qwave(indexQ(:,i)) = ecg4(indexQ(:,i));
end
figure('Name','Q waves','NumberTitle','off'); 
plot(t,Qwave);

TempQ = zeros(N_Qw,1);
for i=1:9
    TempQ = TempQ + Qwave(indexQ(:,i)) ;
end
TempQ = TempQ/9;
figure('Name','Template Q','NumberTitle','off'); 
plot(t(1:N_Qw),TempQ);

stdQ = 0;
noise = stdQ*randn(N,1);
ecg_noise = ecg4 + noise;
[rq,lagsq] = xcorr(ecg_noise,TempQ);
% figure('Name','ECG noise Q','NumberTitle','off'); 
% plot(t,ecg_noise)
% figure('Name','Cross correlation Q','NumberTitle','off'); 
% plot(lagsq/Fs,rq);
% 
% V_thresholdq=0.15;
% hold on
% plot(lagsq/Fs,V_thresholdq*ones(2*N-1,1))
% 
% p_pos = zeros(N,1);

 