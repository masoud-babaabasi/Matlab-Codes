clear
filename = 'data_31_tir.xlsx';
%% single transmiter 1
start = 2;
angle = 0:10:90;
position = -6.5:6.5:26;
power_single_1 = zeros(length(position),length(angle));
efficiency_single_1 = zeros(length(position),length(angle));
for i=1:length(position)
    rangex = 'H'+string(start+ ( i - 1 ) * 10 )+':H'+string(start + i * 10 );
    power_single_1(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'J'+string(start+ ( i - 1 ) * 10 )+':J'+string(start + i * 10 );
    efficiency_single_1(i,:) = table2array(readtable(filename,'Range',rangex));
end
%% plot
figure('Name','single transmiter recieved power')
surf(angle,position,power_single_1)
xlabel('Receiver angle(\circ)')
xlim([0 90])
ylabel('Receiver Position')
zlabel('Received Power(w)')
figure('Name','single transmiter efficiendy')
surf(angle,position,efficiency_single_1)
xlabel('Receiver angle(\circ)')
xlim([0 90])
ylabel('Receiver Position')
zlabel('Efficiency(%)')
%% INphase transmiter 2&4
power_inphase_2_4 = zeros(length(position),length(angle));
efficiency_inphase_2_4 = zeros(length(position),length(angle));
for i=1:length(position)
    rangex = 'U'+string(start+ ( i - 1 ) * 10 )+':U'+string(start + i * 10 );
    power_inphase_2_4(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'W'+string(start+ ( i - 1 ) * 10 )+':W'+string(start + i * 10 );
    efficiency_inphase_2_4(i,:) = table2array(readtable(filename,'Range',rangex));
end
%% plot
figure('Name','single transmiter recieved power')
xlabel('Receiver angle(\circ)')
xlim([0 90])
ylabel('Receiver Position')
zlabel('Received Power(w)')
surf(angle,position,power_inphase_2_4)
figure('Name','single transmiter efficiendy')
surf(angle,position,efficiency_inphase_2_4)
xlabel('Receiver angle(\circ)')
xlim([0 90])
ylabel('Receiver Position')
zlabel('Efficiency(%)')
%% INPHASE 1&2&4
power_dual_inp = zeros(length(position),length(angle));
efficiency__dual_inp = zeros(length(position),length(angle));
for i=1:length(position)
    rangex = 'AJ'+string(start+ ( i - 1 ) * 10 )+':AJ'+string(start + i * 10 );
    power_dual_inp(i,:) = table2array(readtable(filename,'Range',rangex));
        rangex = 'AL'+string(start+ ( i - 1 ) * 10 )+':AL'+string(start + i * 10 );
    efficiency__dual_inp(i,:) = table2array(readtable(filename,'Range',rangex));
end
%% plot
figure('Name','dual inphase recieved power')
surf(angle,position,power_dual_inp)
xlabel('Receiver angle(\circ)')
xlim([0 90])
ylabel('Receiver Position')
zlabel('Received Power(w)')
figure('Name','dual inphase efficiendy')
surf(angle,position,efficiency__dual_inp)
xlabel('Receiver angle(\circ)')
xlim([0 90])
ylabel('Receiver Position')
zlabel('Efficiency(%)')
%% 2&4 INPHASE 1 ANTIPHASE 
power_dual_anti = zeros(length(position),length(angle));
efficiency__dual_anti = zeros(length(position),length(angle));
for i=1:length(position)
    rangex = 'AY'+string(start+ ( i - 1 ) * 10 )+':AY'+string(start + i * 10 );
    power_dual_anti(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'BA'+string(start+ ( i - 1 ) * 10 )+':BA'+string(start + i * 10 );
    efficiency__dual_anti(i,:) = table2array(readtable(filename,'Range',rangex));
end
%% plot 3d
figure('Name','dual antiphase recieved power')
surf(angle,position,power_dual_anti)
xlabel('Receiver angle(\circ)')
xlim([0 90])
ylabel('Receiver Position')
zlabel('Received Power(w)')
figure('Name','dual antiphase efficiendy')
surf(angle,position,efficiency__dual_anti)
xlabel('Receiver angle(\circ)')
xlim([0 90])
ylabel('Receiver Position')
zlabel('Efficiency(%)')
%% plot 2d position angle 0
alpha = 0 ;
figure('Name','power vs position angle =0')
plot( position , power_single_1(: , angle== alpha)*1000 ,'-*');
hold on 
plot( position , power_inphase_2_4(: , angle== alpha)*1000 ,'-o');
plot( position , power_dual_inp(: , angle== alpha)*1000 ,'-square');
xticks(position)
xlim([position(1) position(end)])
xline(0);
ylabel('Received Power(mW)')
legend(["Transmitter 1" "Transmitter 2&3 In-phase" "ALL In-phase"],'Location','northeastoutside')
grid on

figure('Name','effi vs position angle =0')
plot( position , efficiency_single_1(: , angle== alpha) ,'-*');
hold on 
plot( position , efficiency_inphase_2_4(: , angle== alpha) ,'-o');
plot( position , efficiency__dual_inp(: , angle== alpha) ,'-square');
xticks(position)
xlim([position(1) position(end)])
xline(0);
ylabel('Efficiency(%)')
legend(["Transmitter 1" "Transmitter 2&3 In-phase" "ALL In-phase"],'Location','northeastoutside')
grid on

%% plot 2d position angle 90
alpha = 90 ;
figure('Name','power vs position angle =0')
plot( position , power_single_1(: , angle== alpha)*1000 ,'-*');
hold on 
plot( position , power_dual_anti(: , angle== alpha)*1000 ,'-square');
xticks(position)
xlim([position(1) position(end)])
xline(0);
legend(["Transmitter 1" "T2&3 In-phase T1 Anti_phase"],'Location','northeastoutside')
grid on

figure('Name','effi vs position angle =0')
plot( position , efficiency_single_1(: , angle== alpha) ,'-*');
hold on 
plot( position , efficiency__dual_anti(: , angle== alpha) ,'-square');
xticks(position)
xlim([position(1) position(end)])
xline(0);
legend(["Transmitter 1"  "T2&3 In-phase T1 Anti_phase"],'Location','northeastoutside')
grid on
%%
figure('Name','power vs angle position =0')
p = 0 ;
plot( angle , power_dual_inp(position == p , :)*1000 ,'-*');
hold on 
plot( angle , power_dual_anti(position == p , :)*1000 ,'-o');
xlabel('Receiver angle(\circ)')
ylabel('Received Power(mW)')
xlim([0 90])
xticks(angle)
legend(["In-phase" "Anti-phase"],'Location','northeastoutside')
grid on


figure('Name','effi vs angle position =0')
p = 0 ;
plot( angle , efficiency__dual_inp(position == p , :) ,'-*');
hold on 
plot( angle , efficiency__dual_anti(position == p , :) ,'-o');
plot( angle , efficiency_single_1(position == p , :) ,'-s');
plot( angle , efficiency_inphase_2_4(position == p , :) ,'-d');
xlabel('Receiver angle(\circ)')
ylabel('Efficiency(%)')
xticks(angle)
xlim([0 90])
legend(["In-phase" "Anti-phase" "Single" "Dual"],'Location','northeastoutside')
grid on