clear
INPHASE = 1 ;
OFF = 0 ;
ANTIPHASE = 15;
bluetooth_Name = "reciever coil";
s = serialport('COM5' ,115200);
voltage = 2;
%------------------------------
for i= 1:7
    set_transmiter(s , i , OFF);
    pause(.1)
end
f_cpu = 180 ; % MHz
fres = 841 ; % KHz
set_freq(s , f_cpu , fres);
STARTCONVERSIONCMD = zeros(1,2);
STARTCONVERSIONCMD(1) = double('S');
STARTCONVERSIONCMD(2) = double('C');
%%
% set_transmiter(s , transmiter_num , ANTIPHASE);
set_voltage(s, voltage);
voltage = 2:0.2:12;
It = zeros(7,length(voltage));
n_avg = 20;
it_avg = 0;
for transmiter_num = 1:7
    set_transmiter(s , transmiter_num , ANTIPHASE);
    i = 0;
    for v = voltage
        set_voltage(s, v);
%         pause(0.1)
        for ii=1:n_avg
        write(s , STARTCONVERSIONCMD , "uint8");
            read_data = read(s,14,"uint8");
            read_data_16 = read_data(1:2:13) + read_data(2:2:14)*256;
            it_avg = it_avg + read_data_16(transmiter_num);
        end
        i = i + 1;
        It(transmiter_num,i) =it_avg / n_avg;
        it_avg = 0;
    end
    set_transmiter(s , transmiter_num  , OFF);
end
set_voltage(s, 2);
figure('Name','Current Vs voltage');
plot(voltage,It);
%%
set_voltage(s, voltage);
set_freq(s , f_cpu , 857);
voltage = 2:0.2:12;
It_in = zeros(7,7,length(voltage),2);
for transmiter_num1 = 1:7
    set_transmiter(s , transmiter_num1 , ANTIPHASE);
    for transmiter_num2 = 1:7
        if transmiter_num2 == transmiter_num1
            continue;
        end
        set_transmiter(s , transmiter_num2 , ANTIPHASE);
        i = 0;
        for v = voltage
            set_voltage(s, v);
            pause(0.1)
            write(s , STARTCONVERSIONCMD , "uint8");
            read_data = read(s,14,"uint8");
            read_data_16 = read_data(1:2:13) + read_data(2:2:14)*256;
            i = i + 1;
            It_in(transmiter_num1,transmiter_num2,i,:) = read_data_16([transmiter_num1 , transmiter_num2]);
        end
        set_transmiter(s , transmiter_num2  , OFF);
    end
    set_transmiter(s , transmiter_num1 , OFF);
end
set_voltage(s, 2);
%%
%fres dual inphase 857KHz
set_voltage(s, voltage);
for i= 1:7
    set_transmiter(s , i , OFF);
    pause(.1)
end
freq = f_cpu / 2  ./ [100:112] * 1000;
set_transmiter(s , 1 , ANTIPHASE);
pause(0.1)
set_transmiter(s , 7 , ANTIPHASE);
i = 0;
It_dual_in=zeros(length(freq),7);
n_avg = 30;
for fres = freq
    set_freq(s , f_cpu , fres); 
    pause(0.2)
    i = i + 1;
    for j = 1:n_avg
        write(s , STARTCONVERSIONCMD , "uint8");
        read_data = read(s,14,"uint8");
        read_data_16 = read_data(1:2:13) + read_data(2:2:14)*256;
        It_dual_in(i,:) = It_dual_in(i,:) + read_data_16;
    end
    It_dual_in(i,:) = It_dual_in(i,:) / n_avg;
end
figure('Name','Current Vs frequency');
plot(freq,It_dual_in);