clear 

load mutual_k_paper L1 L2 D K R1 
X = D(D/R1 >= -2 & D/R1<=4)/R1;
K1 = K(1,D/R1 >= -2 & D/R1<=4)';
K2 = K(1,D/R1 >= -4 & D/R1<=2)';
% plot(X,K1);
% hold on 
% plot(X,K2);
Rs = 1;
RL = 100;
w = 841e3;
Vs = 1;
M1 = K1 * sqrt(L1*L2);
M2 = K2 * sqrt(L1*L2);
%% single
Z = zeros(length(M1) , 2 , 2);
Z(:,1,:) = [Rs*ones(size(M1)) ,  M1*w*1i];
Z(:,2,:) = [M1*w*1i , -L2*w*1i/RL*ones(size(M1))] ;
I1 = zeros(1,length(M1));
V3 = zeros(1,length(M1));
V = [ Vs ; 0];
for i=1:length(M1)
    I = squeeze(Z(i,:,:))^-1 * V;
    I1(i) = I(1);
    V3(i) = I(2);
end
Pr = abs(V3).^2/RL ; 
Pt = abs(I1) * Vs;
effi1 = Pr ./ Pt;

Z(:,1,:) = [Rs*ones(size(M2)) ,  M2*w*1i];
Z(:,2,:) = [M2*w*1i , -L2*w*1i/RL*ones(size(M2))] ;
I2 = zeros(1,length(M1));
V3 = zeros(1,length(M1));
V = [ Vs ; 0];
for i=1:length(M1)
    I = squeeze(Z(i,:,:))^-1 * V;
    I2(i) = I(1);
    V3(i) = I(2);
end
Pr = abs(V3).^2/RL ; 
Pt = abs(I2) * Vs;
effi2 = Pr ./ Pt;
%% inphase calculation
Z = zeros(length(M1) , 3 , 3);
Z(:,1,:) = [Rs*ones(size(M1)) , zeros(size(M1)) , M1*w*1i];
Z(:,2,:) = [zeros(size(M1)) , Rs*ones(size(M1)) ,  M2*w*1i] ;
Z(:,3,:) = [M1*w*1i ,  M2*w*1i , -L2*w*1i/RL*ones(size(M1))] ;
V = [ Vs ; Vs ; 0];
I1 = zeros(1,length(M1));
I2 = zeros(1,length(M1));
V3 = zeros(1,length(M1));
for i=1:length(M1)
    I = squeeze(Z(i,:,:))^-1 * V;
    I1(i) = I(1);
    I2(i) = I(2);
    V3(i) = I(3);
end
Pr = abs(V3).^2/RL ; 
Pt = (abs(I1) +abs(I2)) * Vs;
effi_inphase = Pr ./ Pt;

figure
plot(X,effi1 * 100);
hold on
plot(X,effi2 * 100);
plot(X,effi_inphase * 100);

xlabel('Normalized position to transmitter Radius')
ylabel('Efficiency(%)')
legend(["Coil1" "Coil2" "In-phase"])
grid on
%%

K1 = K(end,D/R1 >= -2 & D/R1<=4)';
K2 = K(end,D/R1 >= -4 & D/R1<=2)';
M1 = K1 * sqrt(L1*L2);
M2 = K2 * sqrt(L1*L2);
%% anti-phase
Z = zeros(length(M1) , 3 , 3);
Z(:,1,:) = [Rs*ones(size(M1)) , zeros(size(M1)) , M1*w*1i];
Z(:,2,:) = [zeros(size(M1)) , Rs*ones(size(M1)) ,  M2*w*1i] ;
Z(:,3,:) = [M1*w*1i ,  M2*w*1i , -L2*w*1i/RL*ones(size(M1))] ;
I1 = zeros(1,length(M1));
I2 = zeros(1,length(M1));
V3 = zeros(1,length(M1));
V = [ Vs ; Vs ; 0];
for i=1:length(M1)
    I = squeeze(Z(i,:,:))^-1 * V;
    I1(i) = I(1);
    I2(i) = I(2);
    V3(i) = I(3);
end
Pr = abs(V3).^2/RL ; 
Pt = (abs(I1) +abs(I2)) * Vs;
effi_anti = Pr ./ Pt;

figure
plot(X,effi_anti * 100);
%% angle
e = 1e-6;
K1 = K(:,abs(D / R1 - 1) <= e );
K2 = K(:,abs(D / R1  + 1) <= e );
Range_tilt = 0:10:90;
M1 = K1 * sqrt(L1*L2);
M2 = K2 * sqrt(L1*L2);
Z = zeros(length(M1) , 3 , 3);
Z(:,1,:) = [Rs*ones(size(M1)) , zeros(size(M1)) , M1*w*1i];
Z(:,2,:) = [zeros(size(M1)) , Rs*ones(size(M1)) ,  M2*w*1i] ;
Z(:,3,:) = [M1*w*1i ,  M2*w*1i , -L2*w*1i/RL*ones(size(M1))] ;
V = [ Vs ; Vs ; 0];
I1 = zeros(1,length(M1));
I2 = zeros(1,length(M1));
V3 = zeros(1,length(M1));
for i=1:length(M1)
    I = squeeze(Z(i,:,:))^-1 * V;
    I1(i) = I(1);
    I2(i) = I(2);
    V3(i) = I(3);
end
Pr = abs(V3).^2/RL ; 
Pt = (abs(I1) +abs(I2)) * Vs;
effi_inphase_anlge = Pr ./ Pt;

V = [ Vs ; -Vs ; 0];
I1 = zeros(1,length(M1));
I2 = zeros(1,length(M1));
V3 = zeros(1,length(M1));
for i=1:length(M1)
    I = squeeze(Z(i,:,:))^-1 * V;
    I1(i) = I(1);
    I2(i) = I(2);
    V3(i) = I(3);
end
Pr = abs(V3).^2/RL ; 
Pt = (abs(I1) +abs(I2)) * Vs;
effi_anti_anlge = Pr ./ Pt;
figure
plot(Range_tilt , effi_inphase_anlge * 100);
hold on 
plot(Range_tilt , effi_anti_anlge * 100);

xlabel('Tilt angle( \circ )')
ylabel('Efficiency(%)')
legend(["In-phase" "Anti-phase"])
grid on