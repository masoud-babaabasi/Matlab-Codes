clear
L=4.3e-6;
C=3.35e-9;
r = 1;
K = -0.033;
m = K*L;
f_cpu = 180;
disp('L = '+string(L*1e6)+'uH , C = '+string(C*1e9)+'nF , K = '+string(K)+'%')
disp('Fres single = '+string(1/(2*pi*sqrt(L*C))*1e-6)+'MHz'+  ' , microcontroller freq = '+string(f_cpu/round(f_cpu*1e6*(2*pi*sqrt(L*C))/2)/2));
disp('Fres dual Inphase = '+string(1/(2*pi*sqrt((L+m)*C))*1e-6)+'MHz'+  ' , microcontroller freq = '+string(f_cpu/round(f_cpu*1e6*(2*pi*sqrt((L+m)*C))/2)/2));
disp('Fres dual Antiphase = '+string(1/(2*pi*sqrt((L-m)*C))*1e-6)+'MHz'+  ' , microcontroller freq = '+string(f_cpu/round(f_cpu*1e6*(2*pi*sqrt((L-m)*C))/2)/2));
disp('Fres triple Inphase = '+string(1/(2*pi*sqrt((L+2*m)*C))*1e-6)+'MHz'+  ' , microcontroller freq = '+string(f_cpu/round(f_cpu*1e6*(2*pi*sqrt((L+2*m)*C))/2)/2));
a6= C^3*(L-m)*(L+2*m)*(L+3*m);
a5=0;
a4 = -( r^2*C^3*(L+3*m) + C^2*(2*L+m)*(L+3*m) + (L-m)*(L+2*m)*C^2 -r^2*C^3*(2*L+m));
a3 = 0;
a2= (3*L+4*m)*C-r^2*C^2;
a1=0;
a0=-1;
p=[ a6 a5 a4 a3 a2 a1 a0];
freq= roots(p)/(2*pi);
for i=1:length(freq)
    if freq(i) == real(freq(i)) && freq(i) > 0
        fres(i) = freq(i);
        disp('Fres triple Antiphase = '+string(fres(i)*1e-6)+'MHz'+  ' , microcontroller freq = '+string(f_cpu/round(f_cpu*1e6/fres(i)/2)/2))
%         break;
    end
end

% L = 4.3uH , C = 3.35nF , K = -0.033%
% Fres single = 1.3261MHz , microcontroller freq = 1.3235
% Fres dual Inphase = 1.3485MHz , microcontroller freq = 1.3433
% Fres dual Antiphase = 1.3047MHz , microcontroller freq = 1.3043
% Fres triple Inphase = 1.3721MHz , microcontroller freq = 1.3636
% Fres triple Antiphase = 1.3035MHz , microcontroller freq = 1.3043