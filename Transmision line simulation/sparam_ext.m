clear 
filename = 'channel2.s4p';
backplane = sparameters(filename);
data = backplane.Parameters;
freq = backplane.Frequencies;
z0 = backplane.Impedance;
% figure('Name','S parameters Absolute');
% n=4;
% for i=1:n
%     for j=1:n
%         subplot(n,n,i+(j-1)*n)
%         plot(freq*1e-9,20*log10(abs(squeeze(data(j,i,:)))))
%         ylabel('S'+string(j)+string(i))
%         xlabel('Frequency(GHz)')
%         grid on
%     end
% end
% figure('Name','S parameters phase');
% n=2;
% for i=1:n
%     for j=1:n
%         subplot(n,n,i+(j-1)*n)
%         plot(freq*1e-9,20*log10(angle(squeeze(data(j,i,:)))))
%         ylabel('S'+string(j)+string(i))
%         xlabel('Frequency(GHz)')
%         grid on
%     end
% end
%%
Ts = 1e-12;
figure('Name','impulse response');
n = length(freq);
fstp = (freq(end)-freq(1))/(n-1);
fstart = freq(1);
% nfft = 2^(floor(log2(n))+1)*2;
% nfft = 2^16;
fs = 1 / Ts;
nfft = 2/fstp*(fs/2+1);
nfft = max(2^(floor(log2(nfft))+1)*2,2^(floor(log2(n))+1)*2);
f_int = [0:fstp:(fstart/fstp-1)*fstp freq'];
H1=squeeze(data(2,1,:))';
H_int=interp1(freq,H1,f_int(2:end),'linear','extrap');
% H=[1 ones(fstart/fstp-1,1)' conj(squeeze(data(2,1,:)))' zeros(nfft-2*n-2*fstart/fstp+1,1)' squeeze(data(2,1,end:-1:1))' ones(fstart/fstp-1,1)'];
H=[1 conj(H_int) zeros(1,nfft-2*n-2*fstart/fstp+1) H_int(end:-1:1) ];
imp=ifft(H)*length(H);
fs = (fstp * nfft/2-1)*2;
Ts = 1/fs;
t = 0:Ts:(nfft-1)*Ts;
plot(t*1e9,imp)
clear H1 H_int f_int 
%%
Baud = 8e9;
index_t0 = find(imp==max(imp));
ns_baud = floor(fs/Baud);
n_pro = 10;
n_pre = 4;
n_pre = min(n_pre,floor(index_t0/ns_baud));
sbr_sample =[index_t0-(n_pre*ns_baud):ns_baud:index_t0-ns_baud,index_t0:ns_baud:index_t0+n_pro*ns_baud];
t_sbr = t(sbr_sample);
sbr = imp(sbr_sample);
figure('Name','single bit response');
stem(t_sbr*1e9,sbr)
%%
input = dec2bin(0:2^length(sbr)-1)-'0';
input(:,n_pre+1) = 1;
input = 2*input - 1;
for i=1:length(input)
    ISI1(i) = sum(input(i,:).*sbr);
end
input(:,n_pre+1) = -1;
for i=1:length(input)
    ISI0(i) = sum(input(i,:).*sbr);
end
figure('Name','data 1 PDF');
histogram (ISI1,100); 
hold on 
histogram (ISI0,100); 
figure('Name','data 1 CDF');
ecdf(ISI1)
hold on 
[cdf0,x] = ecdf(ISI0);
plot(x,cdf0(end:-1:1))
