%Answer of the second computer assignment of HSSL 2021
%author: Masoud Babaabasi
clear 
filename = 'channel1.s4p';
backplane = sparameters(filename);
% data = backplane.Parameters;
% 
% temp = squeeze(data(2,:,:));
% data(2,:,:) = squeeze(data(3,:,:));
% data(3,:,:) = temp;
% 
% temp = data(:,2,:);
% data(:,2,:) = squeeze(data(:,3,:));
% data(:,3,:) = temp;
% backplane.Parameters = data;
% backplane= cascadesparams(backplane,backplane,backplane);

data = backplane.Parameters;
freq = backplane.Frequencies;
Ts = 1e-12;

H1=squeeze(data(1,2,:));
[imp,t,Ts] = impulse_response(H1,freq,Ts);
% im = conv(imp,imp);
% im = conv(imp,im);
% im = im(1:length(t));
% plot(t,im)
% plot(t*1e9,im)
% imp = im;
% imp = imp/max(imp);
clear im
xlim([0 20])
%% step response
step = cumsum(imp);
step = step / max(step);
figure('Name','step response')
plot(t*1e9,step);
xlim([0 20])
title('step response')
%% Single bit response
Baud = 7e9;
index_t0 = find(imp==max(imp));
ns_baud = floor(1/Ts/Baud);
n_pro = 6;
n_pre = 3;
n_pre = min(n_pre,floor(index_t0/ns_baud));
sbr_sample_index =[index_t0-(n_pre*ns_baud):ns_baud:index_t0-ns_baud  index_t0:ns_baud:index_t0+n_pro*ns_baud];
t_sbr = t(sbr_sample_index);
sbr = imp(sbr_sample_index);
figure('Name','single bit response');
% stem(t_sbr*1e9,sbr)
stem([-n_pre:n_pro],sbr)
ISI = 2*sbr(n_pre+1)-sum(abs(sbr));
%% worst eye
eye_opening = 2*(sbr(n_pre+1) -(sum(abs(sbr))-abs(sbr(n_pre+1))));
disp('eye opening = ' + string(eye_opening));
worst_data = -sign(abs(sbr));
worst_data(n_pre+1) = 1;
figure('Name','wrost data input ')
stem([-n_pre:n_pro],worst_data);

tx_data = kron(worst_data,ones(1,ns_baud));
index = find(t<20e-9+Ts/2 & t>20e-9-Ts/2);
imp = imp(1:index);
t = t(1:index);
rx_data=conv(tx_data,imp);
rx_data = rx_data(1:end-length(tx_data)+1);
rx_data = rx_data/max(abs(rx_data));
figure('Name','wrost data output')
plot(t*1e9,rx_data)
hold on 
plot(t*1e9,-rx_data)
% plot(t(1:length(tx_data))*1e9,tx_data)
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
% histogram (ISI0,100); 
figure('Name','data 1 CDF');
ecdf(ISI1)
hold on 
[cdf0,x] = ecdf(ISI0);
plot(x,cdf0(end:-1:1))
%% FIR transmitter
n_FIR = 5;
l = n_FIR + length(sbr) -1;
H = zeros(l,n_FIR);
for i = 1:l
    if i >= n_FIR && i <= l - n_FIR + 1
    H(i,:) = sbr(i:-1:i-n_FIR+1);
    elseif i < n_FIR
    H(i,1:i) = sbr(i:-1:1);
    else 
    H(i,i-l+n_FIR:n_FIR) = sbr(length(sbr):-1:i-n_FIR+1);   
    end
end
Ydes = zeros(1,l);
Ydes((n_pre+1)+1) = 1;
Wls =((inv(H'*H)*H')*Ydes')';
Wls = Wls / norm(Wls);
Y = H * Wls';
figure('Name','single bit response')
subplot(2,1,1)
stem(sbr);
subplot(2,1,2)
stem(Y);

Wls_conv = kron(Wls,ones(1,ns_baud));
imp_eq = conv(imp,Wls_conv);
rx_data_eqe=conv(tx_data,imp_eq);
rx_data_eqe = rx_data_eqe(1:end-length(tx_data)+1);
rx_data_eqe = rx_data_eqe/max(abs(rx_data_eqe));
figure('Name','wrost data output')
plot(rx_data_eqe)
hold on 
% plot(t(1:length(tx_data))*1e9,tx_data)
plot(-rx_data_eqe)

index_t0 = find(imp_eq==max(imp_eq));
n_pre = min(n_pre,floor(index_t0/ns_baud));
sbr_sample_index =[index_t0-(n_pre*ns_baud):ns_baud:index_t0-ns_baud  index_t0:ns_baud:index_t0+n_pro*ns_baud];
sbr_equ = imp_eq(sbr_sample_index);
sbr_equ = sbr_equ/abs(max(sbr_equ));
eye_opening_equ = 2*(sbr_equ(n_pre+1) -(sum(abs(sbr_equ))-abs(sbr_equ(n_pre+1))));
disp('eye opening equalized = ' + string(eye_opening_equ));

%% eye
n_bits = 1000;
bit_stream = 2*round(rand(1,n_bits))-1;
random_tx = kron(bit_stream,ones(1,ns_baud));
random_rx = conv(random_tx,imp);
random_rx = random_rx / max(abs(random_rx));
random_rx_eq = conv(random_tx,imp_eq);
random_rx_eq = random_rx_eq / max(abs(random_rx_eq));
index_t0 = find(imp==max(imp));
index_t0_eq = find(imp_eq==max(imp_eq));
t_delay = t(index_t0);
t_delay_eq = t(index_t0_eq);
figure('Name','eye diagram')
r2 = round(t_delay/Ts)+1+120;
r2_eq = round(t_delay_eq/Ts)+1;
t_unit = linspace(0,1/Baud,ns_baud)*1e12;
for i=0:n_bits-1
    subplot(2,1,1)
%     r1 = round((t_delay+1/Baud*i)/Ts);
%     r2 = round((t_delay+1/Baud*(i+1))/Ts);
   r1 = r2+1;
   r2 = r1+ns_baud-1;
   plot(t_unit, random_rx(r1:r2))
   hold on
   subplot(2,1,2)
%     r1 = round((t_delay_eq+1/Baud*i)/Ts)+1;
%     r2 = round((t_delay_eq+1/Baud*(i+1))/Ts)+1;
   r1_eq = r2_eq+1; 
   r2_eq = r1_eq +ns_baud-1;
   plot(t_unit, random_rx_eq(r1_eq:r2_eq))
   hold on
end
subplot(2,1,1)
title(string(Baud*1e-9)+'Gbps chennel eye diagram')
xlabel('ps')
xlim([0 t_unit(end)])
subplot(2,1,2)
title(string(Baud*1e-9)+'Gbps FIR equalized chennel eye diagram')
xlabel('ps')
xlim([0 t_unit(end)])
%% DEF
n_dfe = 3;
din_dfe = zeros(1,n_dfe);
din_dfe_eq = zeros(1,n_dfe);
rx_dfe = -ones(1,ns_baud);
rx_dfe_eq = -ones(1,ns_baud);
figure('Name','eye diagram dfe')
r2 = round(t_delay/Ts)+1;
r2_eq = round(t_delay_eq/Ts)+1;
sbr_equ = sbr_equ / sbr_equ(n_pre+1);
for i =0:n_bits-1
%     r1 = round((t_delay+1/Baud*i)/Ts)+1;
%     r2 = round((t_delay+1/Baud*(i+1))/Ts)+1;
    subplot(2,1,1)
    r1 = r2+1;
    r2 = r1+ns_baud-1;
    rx_dfe = random_rx(r1:r2)-sum(din_dfe(1:end).*sbr(n_pre+2:n_pre+n_dfe+1));
    din_dfe(2:end) = din_dfe(1:end-1);
    din_dfe(1) = sign(rx_dfe(1));    
    plot(rx_dfe);
    hold on
    
    subplot(2,1,2)
    r1_eq = r2_eq+1;
    r2_eq = r1_eq+ns_baud-1;
    rx_dfe_eq = random_rx_eq(r1_eq:r2_eq)-sum(din_dfe_eq(1:end).*sbr_equ(n_pre+2:n_pre+n_dfe+1));
    din_dfe_eq(2:end) = din_dfe_eq(1:end-1);
    din_dfe_eq(1) = sign(rx_dfe_eq(1));    
    plot(rx_dfe_eq);
    hold on
end
subplot(2,1,1)
title(string(Baud*1e-9)+'Gbps chennel eye diagram')
xlabel('ps')
xlim([0 t_unit(end)])
subplot(2,1,2)
title(string(Baud*1e-9)+'Gbps DFE equalized chennel eye diagram')
xlabel('ps')
xlim([0 t_unit(end)])
%% 4 level eye
bit_stream2 = 2*bit_stream(2:2:end) + bit_stream(1:2:end-1); 
random_tx2 = kron(bit_stream2,ones(1,ns_baud));
random_rx2 = conv(random_tx2,imp);
random_rx2 = random_rx2 / max(abs(random_rx2));
random_rx_eq2 = conv(random_tx2,imp_eq);
random_rx_eq2 = random_rx_eq2 / max(abs(random_rx_eq2));
figure('Name','4 level eye diagram')
r2 = round(t_delay/Ts)+1;
r2_eq = round(t_delay_eq/Ts)+1;
for i=0:n_bits/2-1
    subplot(2,1,1)
%     r1 = round((t_delay+1/Baud*i)/Ts);
%     r2 = round((t_delay+1/Baud*(i+1))/Ts);
   r1 = r2+1;
   r2 = r1+ns_baud-1;
   plot(t_unit, random_rx2(r1:r2))
   hold on
   subplot(2,1,2)
%     r1 = round((t_delay_eq+1/Baud*i)/Ts)+1;
%     r2 = round((t_delay_eq+1/Baud*(i+1))/Ts)+1;
   r1_eq = r2_eq+1; 
   r2_eq = r1_eq +ns_baud-1;
   plot(t_unit, random_rx_eq2(r1_eq:r2_eq))
   hold on
end

subplot(2,1,1)
title(string(2*Baud*1e-9)+'Gbps chennel eye diagram')
xlabel('ps')
xlim([0 t_unit(end)])
subplot(2,1,2)
title(string(2*Baud*1e-9)+'Gbps FIR equalized chennel eye diagram')
xlabel('ps')
xlim([0 t_unit(end)])