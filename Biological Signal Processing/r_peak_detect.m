%ECG file load
load('103m.mat');
%enter duration of the ECG signal in seconds
t_sampling=60;
%lead 1 data
s1 = val(1,:);
%lead 2 data
s2 = val(2,:);
%calculation of sampling frequency
fs = length(s1) / t_sampling;

t=0:(1/fs):(t_sampling-(1/fs));
figure('Name','ECG Data','NumberTitle','off');
subplot(2,1,1)
plot(t,s1)
title('lead1')
subplot(2,1,2)
plot(t,s2)
title('lead2')
%*******Heart Rate Variablity calculation
%enter a threshold level that all R peaks are above this threshold
V_threshold=100;
%****************
f=0;
first_peak=1;
k=1;
for i=1:length(s1)
   if( s1(i) > V_threshold &&  f == 0 )
    R1=i;
    f=1;
   end
   if( s1(i) <= V_threshold &&  f == 1 )
       R2=i;
       f=0;
       Vmax = s1(R1);
       I_Vmax2=R1;
       for j=R1+1:R2
           if( s1(j) > Vmax )
               Vmax = s1(j);
                I_Vmax2=j;
           end
       end
       if( first_peak == 1 )
          I_Vmax1 = I_Vmax2;
          first_peak=0;
       else
           hrv(k) = fs * 60 / (I_Vmax2 - I_Vmax1);
           t1(k) = j / fs;
           k = k + 1;
           I_Vmax1 = I_Vmax2;
       end
   end
end
figure('Name','Heart rate variablity','NumberTitle','off');
plot(t1,hrv)
title('HRV')
%*************
% t_realtime=3;
% figure('Name','Real time lead1','NumberTitle','off');
% for i=1:(t_sampling-t_realtime/2)*fs
% plot(t(i:i+t_realtime/2*fs),s1(i:i+t_realtime/2*fs))
% xlim([t(i) t(i+t_realtime*fs)])
% ylim([min(s1) max(s1)])
% pause(1/fs)
% end

% t_realtime=3;
% sig=zeros(1,t_realtime*fs);
% figure('Name','Real time2 lead1','NumberTitle','off');
% for i=1:(t_sampling-t_realtime)
%     for j=1:t_realtime*fs
% %         tic
%         sig(j)=s1((i-1)*t_realtime*fs+j);
%         plot(t((i-1)*t_realtime*fs+1:i*t_realtime*fs),sig)
%         xlim([t((i-1)*t_realtime*fs+1) t(i*t_realtime*fs)])
%         ylim([min(s1) max(s1)])
% %         t_ex=toc;
% %         pause(1/fs-t_ex)
%     end
% end