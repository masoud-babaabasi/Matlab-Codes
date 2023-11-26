%Answer of the third computer assignment of HSSL 2021
%author: Masoud Babaabasi
clear 
filename = 'channel1.s4p';
backplane = sparameters(filename);
data = backplane.Parameters;
freq = backplane.Frequencies;
Ts = 1e-12;

H1=squeeze(data(1,2,:));
[imp,t,Ts] = impulse_response(H1,freq,Ts);
imp = imp/max(imp);
clear im
xlim([0 10])
%% Single bit response
Baud = 10e9;
index_t0 = find(imp==max(imp));
ns_baud = floor(1/Ts/Baud);
n_pro = 6;
n_pre = 3;
n_pre = min(n_pre,floor(index_t0/ns_baud));
sbr_sample_index =[index_t0-(n_pre*ns_baud):ns_baud:index_t0-ns_baud  index_t0:ns_baud:index_t0+n_pro*ns_baud];
sbr = imp(sbr_sample_index);
ISI = 2*sbr(n_pre+1)-sum(abs(sbr));
eye_opening = 2*(sbr(n_pre+1) -(sum(abs(sbr))-abs(sbr(n_pre+1))));
disp('eye opening at ' + string(Baud*1e-9) +'GHz = '+ string(eye_opening));
%% FIR transmitter
n_FIR = 4;
n_pre_FIR = 1;
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
Ydes((n_pre+1)+n_pre_FIR) = 1;
Wls =((inv(H'*H)*H')*Ydes')';
% Wls = Wls / norm(Wls);
Wls = Wls / sum(abs(Wls));

Wls_conv = kron(Wls,ones(1,ns_baud));
imp_eq = conv(imp,Wls_conv);
index_t0 = find(imp_eq==max(imp_eq));
n_pre = min(n_pre,floor(index_t0/ns_baud));
sbr_sample_index =[index_t0-(n_pre*ns_baud):ns_baud:index_t0-ns_baud  index_t0:ns_baud:index_t0+n_pro*ns_baud];
sbr_equ = imp_eq(sbr_sample_index);
sbr_equ = sbr_equ/abs(max(sbr_equ));
eye_opening_equ = 2*(sbr_equ(n_pre+1) -(sum(abs(sbr_equ))-abs(sbr_equ(n_pre+1))));
disp('eye opening equalized at '+string(Baud*1e-9)+'GHz = ' + string(eye_opening_equ));
figure('Name','equalized SBR')
subplot(2,1,1)
stem(sbr)
title('original SBR')
subplot(2,1,2)
stem(sbr_equ)
title('FIR sbr')
disp('FIR at '+string(Baud*1e-9)+' GHz = '+strjoin(string(Wls),' , '))