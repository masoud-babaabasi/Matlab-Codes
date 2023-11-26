N = 15; FS = 400; fc = 50; Nfft = 2^16;
h1 = fir2(N-1,[0 50/(FS/2) 50/(FS/2) 100/(FS/2) 200/(FS/2)],[1 1  0.01 0.005 0] ,hamming(N));
h2 = fir2(N-1,[0 50/(FS/2) 50/(FS/2) 100/(FS/2) 200/(FS/2)],[1 1  0.01 0.005 0] ,rectwin(N));
f = FS/Nfft:FS/Nfft:FS/2;
Hf1 = 20*log10(abs(fft(h1,Nfft)));
Hf2 = 20*log10(abs(fft(h2,Nfft)));
figure('Name','frequency domain H(f)','NumberTitle','off');
subplot(2,1,1);
plot(f,Hf1(1:Nfft/2));
set(gca,'XTick',(0:25:200))
set(gca,'YTick',(-100:20:0))
ylim([-100 0]);
title('Hamming window');
grid on;
subplot(2,1,2);
plot(f,Hf2(1:Nfft/2));
title('Rectangle window');
ylim([min(Hf2) max(Hf2)]);
set(gca,'XTick',(0:25:200))
set(gca,'YTick',(-100:20:0))
ylim([-100 0]);
grid on;
