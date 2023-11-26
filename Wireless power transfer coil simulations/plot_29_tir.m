clear
filename = 'data_29_tir.xlsx';
%% single transmiter 4
start = 2;
angle = 0:10:90;
position = [0:6:48];
power_single_4 = zeros(length(position),length(angle));
efficiency_single_4 = zeros(length(position),length(angle));
for i=1:length(position)
    rangex = 'H'+string(start+ ( i - 1 ) * 10 )+':H'+string(start + i * 10 );
    power_single_4(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'J'+string(start+ ( i - 1 ) * 10 )+':J'+string(start + i * 10 );
    efficiency_single_4(i,:) = table2array(readtable(filename,'Range',rangex));
end
%% plot 3d
figure('Name','single transmiter recieved power')
surf(angle,position,power_single_4*1000)
xlabel('Receiver angle(\circ)')
xlim([0 90])
ylabel('Receiver Position')
zlabel('Received Power(mW)')
figure('Name','single transmiter efficiendy')
surf(angle,position,efficiency_single_4)
xlabel('Receiver angle(\circ)')
xlim([0 90])
ylabel('Receiver Position')
zlabel('Efficiency(%)')
%% single transmiter 2
position = [0:6:48]-24;
power_single_2 = zeros(length(position),length(angle));
efficiency_single_2 = zeros(length(position),length(angle));
for i=1:length(position)
    rangex = 'S'+string(start+ ( i - 1 ) * 10 )+':S'+string(start + i * 10 );
    power_single_2(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'U'+string(start+ ( i - 1 ) * 10 )+':U'+string(start + i * 10 );
    efficiency_single_2(i,:) = table2array(readtable(filename,'Range',rangex));
end
%% plot 3d
figure('Name','single transmiter recieved power')
xlabel('Receiver angle(\circ)')
xlim([0 90])
ylabel('Receiver Position')
zlabel('Received Power(w)')
surf(angle,position,power_single_2)
figure('Name','single transmiter efficiendy')
surf(angle,position,efficiency_single_2)
xlabel('Receiver angle(\circ)')
xlim([0 90])
ylabel('Receiver Position')
zlabel('Efficiency(%)')
%% INPHASE
power_dual_inp = zeros(length(position),length(angle));
efficiency__dual_inp = zeros(length(position),length(angle));
for i=1:length(position)
    rangex = 'AF'+string(start+ ( i - 1 ) * 10 )+':AF'+string(start + i * 10 );
    power_dual_inp(i,:) = table2array(readtable(filename,'Range',rangex));
        rangex = 'AH'+string(start+ ( i - 1 ) * 10 )+':AH'+string(start + i * 10 );
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
%% ANTIPHASE
power_dual_anti = zeros(length(position),length(angle));
efficiency__dual_anti = zeros(length(position),length(angle));
for i=1:length(position)
    rangex = 'AS'+string(start+ ( i - 1 ) * 10 )+':AS'+string(start + i * 10 );
    power_dual_anti(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'AU'+string(start+ ( i - 1 ) * 10 )+':AU'+string(start + i * 10 );
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

%% plot 2d
figure('Name','power vs position angle =0')
alpha = 0 ;
% subplot(2,1,1)
plot( position , power_single_4(: , angle== alpha)*1000 ,'-square');
hold on 
plot( position , power_single_2(: , angle== alpha)*1000 ,'-o');
plot( position , power_dual_inp(: , angle== alpha)*1000 ,'-*');
xlabel('Receiver Position')
ylabel('Received Power(mW)')
ax = gca;
xticks(position)
xlim([position(1) position(end)])
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
xline(0);
legend(["Transmitter 1" "Transmitter 2" "Both Inphase"],'Location','northeastoutside')
grid on

% subplot(2,1,2)
figure('Name','Efficiency vs position angle =0')
plot( position , efficiency_single_4(: , angle== alpha) ,'-square');
hold on 
plot( position , efficiency_single_2(: , angle== alpha) ,'-o');
plot( position , efficiency__dual_inp(: , angle== alpha) ,'-*');
plot( position , efficiency__dual_anti(: , angle== alpha) ,'-d');
xlabel('Receiver Position(mm)')
ylabel('Efficiency(%)')
xticks(position)
xlim([position(1) position(end)])
xline(0);
legend(["Transmitter 1" "Transmitter 2" "In-phase" "Anti-phase"],'Location','northeastoutside')
grid on

% figure('Name','FOM vs position angle =0')
% plot( position , efficiency_single_4(: , angle== alpha).*power_single_4(: , angle== alpha) ,'-square');
% hold on 
% plot( position , efficiency_single_2(: , angle== alpha).*power_single_2(: , angle== alpha) ,'-o');
% plot( position , efficiency__dual_inp(: , angle== alpha).*power_dual_inp(: , angle== alpha) ,'-*');
% xlabel('Receiver Position')
% ylabel('FOM')
% xticks(position)
% xlim([position(1) position(end)])
% xline(0);
% legend(["Transmitter 1" "Transmitter 2" "Both Inphase"],'Location','northeastoutside')
% grid on
%% plot 2d
% figure('Name','power vs position angle =90')
alpha = 90 ;
% plot( position , power_single_4(: , angle== alpha)*1000 ,'-square');
% hold on 
% plot( position , power_single_2(: , angle== alpha)*1000 ,'-o');
% plot( position , power_dual_anti(: , angle== alpha)*1000 ,'-*');
% xlabel('Receiver Position')
% ylabel('Received Power(mW)')
% xticks(position)
% xlim([position(1) position(end)])
% ax = gca;
% % ax.XAxisLocation = 'origin';
% % ax.YAxisLocation = 'origin';
% xline(0);
% legend(["Transmitter 1" "Transmitter 2" "Anti-phase"],'Location','northeastoutside')
% grid on

figure('Name','Efficiency vs position angle =90')
plot( position , efficiency_single_4(: , angle== alpha) ,'-square');
hold on 
plot( position , efficiency_single_2(: , angle== alpha) ,'-o');
plot( position , efficiency__dual_inp(: , angle== alpha) ,'-*');
plot( position , efficiency__dual_anti(: , angle== alpha) ,'-d');
xlabel('Receiver Position(mm)')
ylabel('Efficiency(%)')
xticks(position)
xlim([position(1) position(end)])
xline(0);
legend(["Transmitter 1" "Transmitter 2" "In-phase" "Anti-phase"],'Location','northeastoutside')
grid on
%%
% figure('Name','power vs angle position =0')
% p = 0 ;
% plot( angle , power_dual_inp(position == p , :)*1000 ,'-*');
% hold on 
% plot( angle , power_dual_anti(position == p , :)*1000 ,'-o');
% plot( angle , power_single_4(position == p , :)*1000 ,'-s');
% plot( angle , power_single_2(position == p , :)*1000 ,'-d');
% xlabel('Receiver angle(\circ)')
% ylabel('Received Power(mW)')
% xlim([0 90])
% xticks(0:10:90);
% legend(["In-phase" "Anti-phase" "Coil4" "Coil2"],'Location','northeastoutside')
% grid on

figure('Name','effi vs angle position =0')
p = 0 ;
plot( angle , efficiency_single_4 (position == p , :) ,'-s');
hold on 
plot( angle , efficiency_single_2 (position == p , :) ,'-o');
plot( angle , efficiency__dual_inp(position == p , :) ,'-*');
plot( angle , efficiency__dual_anti(position == p , :) ,'-d');


xlabel('Receiver angle(\circ)')
ylabel('Efficiency(%)')
xlim([0 90])
xticks(0:10:90);
legend(["Transmitter 1" "Transmitter 2" "In-phase" "Anti-phase"],'Location','northeastoutside')
grid on