clear
INPHASE = 1 ;
OFF = 0 ;
ANTIPHASE = 15;
bluetooth_Name = "reciever coil";
load('calibrate_currents_vs_voltage_noload_avg.mat','It')
load('calibrate_currents_vs_voltage_noload_avg.mat','voltage')
s = serialport('COM5' ,115200);
Vt = 5;
STARTCONVERSIONCMD = zeros(1,2);
STARTCONVERSIONCMD(1) = double('S');
STARTCONVERSIONCMD(2) = double('C');
%------------------------------
f_cpu = 180 ; % MHz
fres = 844 ; % KHz
set_freq(s , f_cpu , fres)
set_voltage(s, Vt);
for i= 1:7
    set_transmiter(s , i , OFF);
    pause(.1)
end
n_avg =5;
i_avg = 0;
for transmiter_num = 1:7
    set_transmiter(s , transmiter_num , ANTIPHASE);
%     pause(0.1);
    for ii=1:n_avg
        write(s , STARTCONVERSIONCMD , "uint8");
        read_data = read(s,14,"uint8");
        read_data_16 = read_data(1:2:13) + read_data(2:2:14)*256;
        i_avg = i_avg + read_data_16(transmiter_num);
    end
    dev(transmiter_num) = (It(transmiter_num, Vt==voltage ) - i_avg/n_avg)/It(transmiter_num, Vt==voltage );
    set_transmiter(s , transmiter_num , OFF);
    i_avg = 0;
end
[m , transmiter_num] = max(dev);
set_transmiter(s , transmiter_num , ANTIPHASE);
blutoth_live = 0;
%%
while(blutoth_live ~= 1)
    pause(0.2);
    list = blelist;
    n_list = size(list,1);
    for i = 1:n_list
        if strcmp(bluetooth_Name,list{i,2}) == 1
            blutoth_live = 1;
            break;
        end
    end
    if blutoth_live ~= 1
        Vt = Vt + 0.5 ;
        if Vt > 12 
            Vt = 12;
        end
        set_voltage(s, Vt);
    end
end
% Vt = 5;
% set_voltage(s , Vt);
%%
bluetooth_Name = "reciever coil";
b = ble(bluetooth_Name);
charac = b.Characteristics;
ServUUID = charac(5,2);
CharUUID = charac(5,4);
c = characteristic(b,ServUUID{:,:},CharUUID{:,:});
ble_data = read(c);
ble_data = read(c);
pause(0.5);
clear b c charac
b = ble(bluetooth_Name);
charac = b.Characteristics;
ServUUID = charac(5,2);
CharUUID = charac(5,4);
c = characteristic(b,ServUUID{:,:},CharUUID{:,:});

% s = serialport('COM5' ,115200);
ts =50*1e-3;
Twindow = 5;
t = 0:ts:Twindow;
n = length(t);
data = zeros(n,11);
i = 1;
y_max = 2200;
figure('Name','Pot values');
 for j=1:11
        scrollsubplot(10,1,[3*j-2,3*j]);
        h(j) = plot(t,data(:,j));
        if(j <= 7 )
            txt(j)=text(0.1,400,['adc = ',num2str(data(n,j))]);
            xlim([0 Twindow]);
            xticks([0:1:Twindow])
            ylim([0 y_max]);
            ylabel("coil"+string(j),'fontsize',14);
            grid on;
            %set(get(gca,'ylabel'),'rotation',0)
            yticks([0:200:y_max]);
        elseif( j== 8 )
           txt(j)=text(0.1,3,['Vr = ',num2str(data(n,j))]);
            xlim([0 Twindow]);
            xticks([0:1:Twindow])
            ylim([0 12]);
            ylabel("V receiver",'fontsize',14);
            grid on;
            %set(get(gca,'ylabel'),'rotation',0)
            yticks([0:2:12]); 
        elseif( j== 9 )
           txt(j)=text(0.1,20,['Ir = ',num2str(data(n,j))]);
            xlim([0 Twindow]);
            xticks([0:1:Twindow])
            ylim([0 70]);
            ylabel("I receiver",'fontsize',14);
            grid on;
            %set(get(gca,'ylabel'),'rotation',0)
            yticks([0:10:70]);
        elseif( j== 10 )
           txt(j)=text(0.1,125,['Pr = ',num2str(data(n,j))]);
            xlim([0 Twindow]);
            xticks([0:1:Twindow])
            ylim([0 500]);
            ylabel("P receiver",'fontsize',14);
            grid on;
            %set(get(gca,'ylabel'),'rotation',0)
            yticks([0:50:500]); 
        elseif( j== 11 )
           txt(j)=text(0.1,2,['V transmiter = ',num2str(data(n,j))]);
            xlim([0 Twindow]);
            xticks([0:1:Twindow])
            ylim([0 12]);
            ylabel("V transmiter",'fontsize',14);
            grid on;
            %set(get(gca,'ylabel'),'rotation',0)
            yticks([0:2:12]); 
        end
 end
drawnow;
%%
STARTCONVERSIONCMD = zeros(1,2);
STARTCONVERSIONCMD(1) = double('S');
STARTCONVERSIONCMD(2) = double('C');
read_data = zeros(1,14);
read_data_16 = zeros(1,7);
time_count = 0;
t_start = tic;
BLE_data_in = zeros(1,2);
t_delay = 0.2;%system delay for feedback 
Vref = 8 ;
dead_zone = 0.1 ; 
Kp = 0.05;
while(1)
    ble_data = read(c);
%     if length(ble_data) ~= 2
%         clear ble_data;
%         continue;
%     end
%     BLE_data_in =sscanf(string(char(ble_data)),"vin=%f , Iin = %fmA");
    BLE_data_in(1) = ble_data(1) * 14 / 255 + 3 ;
    BLE_data_in(2) = ble_data(2) * 60 / 255 + 10 ;
    write(s , STARTCONVERSIONCMD , "uint8");
    read_data = read(s,14,"uint8");
    read_data_16 = read_data(1:2:13) + read_data(2:2:14)*256;
    data(1:n-1,:) = data(2:n,:);
    data(n,1:7) = read_data_16;
    data(n,8) = BLE_data_in(1);
    data(n,9) = BLE_data_in(2);
    data(n,10) = BLE_data_in(1) * BLE_data_in(2);
    data(n,11) = Vt;
    [m , index_v] = min(abs(Vt-voltage));
    deviation_n = (It(transmiter_num, index_v ) - read_data_16(transmiter_num))/It(transmiter_num, index_v ) ;
    if deviation_n < 0.6 * dev(transmiter_num)
        set_transmiter(s , transmiter_num , OFF);
        for transmiter_num = 1:7
            set_transmiter(s , transmiter_num , ANTIPHASE);
        %     pause(0.001);
            write(s , STARTCONVERSIONCMD , "uint8");
            read_data = read(s,14,"uint8");
            read_data_16 = read_data(1:2:13) + read_data(2:2:14)*256;
            dev(transmiter_num) = (It(transmiter_num, index_v ) - read_data_16(transmiter_num))/It(transmiter_num, index_v );
            set_transmiter(s , transmiter_num , OFF);
        end
        [m , transmiter_num] = max(dev);
        set_transmiter(s , transmiter_num , ANTIPHASE);
    end
    for j=1:11
        h(j).YData = data(:,j);
        if( j <= 7 )
            txt(j).String = ['mean adc = ',num2str(round(mean(data(:,j))))];
        elseif( j == 8 )
            t_stop = toc(t_start);
            if t_stop >= t_delay
                time_count = 0;
            end
            vr = data(n,j);
            txt(j).String = ['Vr  = ', num2str(vr)];
%             if(vr <= 9 && time_count == 0)
%                 t_start = tic;
%                 voltage = voltage + 0.1;
%                 set_voltage(s , voltage);
%                 time_count = 1;
%             elseif(vr >= 10 && time_count == 0)
%                 t_start = tic;
%                 voltage = voltage - 0.1;
%                 set_voltage(s , voltage);
%                 time_count = 1;
%             elseif time_count == 0
%                 t_start = tic;
%             end
            if abs( vr - Vref ) >= dead_zone && time_count == 0
                t_start = tic;
                Vt = Vt - Kp * ( vr - Vref );
                if Vt > 12 
                    Vt = 12 ;
                elseif Vt < 2
                    Vt = 2;
                end
                set_voltage(s , Vt);
                time_count = 1;  
            elseif time_count == 0
                t_start = tic;
            end
        elseif( j == 9 )
             txt(j).String = ['Ir = ', num2str(data(n,j)) ];
        elseif( j == 10 )
            txt(j).String = ['P receiver  = ', num2str(data(n,j))];
        elseif( j == 11 )
            txt(j).String = ['V transmiter  = ', num2str(data(n,j))];
        end        
    end
    %pause(ts)
    i = i + 1;
end

