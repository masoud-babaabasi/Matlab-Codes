function [imp,t,Ts] = impulse_response(Hs,freq,Ts_des)
%This code will determine the inverse fourier transform of the given Transfer
%function based on the frequency information given and the desired
%samppling time
%author: Masoud Babaabasi @ March 2021
n = length(freq);
fstp = (freq(end)-freq(1))/(n-1);
fstart = freq(1);
fs = 1 / Ts_des;
nfft = 2/fstp*(fs/2+1);
nfft = max(2^(floor(log2(nfft))+1)*2,2^(floor(log2(n))+1)*2);
size_in = size(freq);
if(size_in(1) ~= 1) 
    freq = freq';
end
f_int = [0:fstp:(fstart/fstp-1)*fstp freq];
size_in = size(Hs);
if(size_in(1) ~= 1) 
    Hs = Hs';
end
H_int=interp1([0 freq],[1 Hs],f_int,'spline');
H_int = H_int(2:end);
H=[1 conj(H_int) zeros(1,nfft-2*n-2*fstart/fstp+1) H_int(end:-1:1) ];
% plot(abs(H))
imp=ifft(H)*length(H);
imp = imp / max(imp);
fs = (fstp * nfft/2-1)*2;
Ts = 1/fs;
t = 0:Ts:(nfft-1)*Ts;
figure('Name','Impulse')
plot(t*1e9,imp)
xlabel('Time (ns)')
title('Normalized impulse response')
disp('calculated number of FFT Points is ' + string(nfft)+'(2^'+string(log2(nfft))+')')
disp('calculated Sampling time due to FFT points (Ts) is ' + string(Ts*1e12)+'ps')
end

