clear
filename = 'data_6&7tir.xlsx';
%%
start = 117;
power = zeros(11,10);
angle = 0:10:90;
position = -25:5:25;
for i=1:11
    rangex = 'D'+string(start+ ( i - 1 ) * 10 )+':D'+string(start + i * 10 );
    power(i,:) = table2array(readtable(filename,'Range',rangex));
end
figure('Name','single transmiter recieved power')
surf(angle,position,power)
%%
filename = 'data_8tir.xlsx';
start = 12;
ADC1_in = zeros(4,10);
ADC2_in = zeros(4,10);
angle = 0:10:90;
position = 10:5:25;
for i=1:4
    rangex = 'N'+string(start+ ( i - 1 ) * 10 )+':N'+string(start + i * 10 );
    ADC1_in(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'O'+string(start+ ( i - 1 ) * 10 )+':O'+string(start + i * 10 );
    ADC2_in(i,:) = table2array(readtable(filename,'Range',rangex));
end
figure('Name','dual transmiter currents')
surf(angle,position,ADC1_in)
hold on
surf(angle,position,ADC2_in)
%%
filename = 'data_9tir.xlsx';
start = 3;
power = zeros(11,10);
curent = zeros(11,10);
angle = 0:10:90;
position = -25:5:25;
for i=1:11
    rangex = 'F'+string(start+ ( i - 1 ) * 10 )+':F'+string(start + i * 10 );
    power(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'C'+string(start+ ( i - 1 ) * 10 )+':C'+string(start + i * 10 );
    curent(i,:) = table2array(readtable(filename,'Range',rangex));
end
figure('Name','single transmiter recieved power 9 tir')
surf(angle,position,power)
legend('Power')
xlabel('angle')
ylabel('position')
figure('Name','single transmiter current 9 tir')
surf(angle,position,curent)
legend('current')
xlabel('angle')
ylabel('position')
%%
power = zeros(11,10);
curent1 = zeros(11,10);
curent2 = zeros(11,10);
angle = 0:10:90;
position = 0:5:50;
for i=1:11
    rangex = 'N'+string(start+ ( i - 1 ) * 10 )+':N'+string(start + i * 10 );
    power(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'J'+string(start+ ( i - 1 ) * 10 )+':J'+string(start + i * 10 );
    curent1(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'K'+string(start+ ( i - 1 ) * 10 )+':K'+string(start + i * 10 );
    curent2(i,:) = table2array(readtable(filename,'Range',rangex));
end
figure('Name','dual transmiter recieved power INPHASE 9 tir')
surf(angle,position,power)
legend('Power INPHASE')
xlabel('angle')
ylabel('position')
figure('Name','dual transmiter currents INPHASE 9 tir')
CO(:,:,3) = zeros(size(power)); % red
CO(:,:,2) = zeros(size(power)); % green
CO(:,:,1) = ones(size(power)); % blue
surf(angle,position,curent1,CO,'FaceColor','interp')
hold on 
CO(:,:,3) = zeros(size(power)); % red
CO(:,:,2) = ones(size(power)); % green
CO(:,:,1) = zeros(size(power)); % blue
surf(angle,position,curent2,CO,'FaceColor','interp')
legend('I1 INPHASE','I2 INPHASE')
xlabel('angle')
ylabel('position')
%%
power = zeros(11,10);
curent1 = zeros(11,10);
curent2 = zeros(11,10);
angle = 0:10:90;
position = 0:5:50;
for i=1:11
    rangex = 'W'+string(start+ ( i - 1 ) * 10 )+':W'+string(start + i * 10 );
    power(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'S'+string(start+ ( i - 1 ) * 10 )+':S'+string(start + i * 10 );
    curent1(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'T'+string(start+ ( i - 1 ) * 10 )+':T'+string(start + i * 10 );
    curent2(i,:) = table2array(readtable(filename,'Range',rangex));
end
figure('Name','dual transmiter recieved power ANTIPHASE 9 tir')
surf(angle,position,power)
legend('Power ANTIPHASE')
xlabel('angle')
ylabel('position')
figure('Name','dual transmiter currents ANTIPHASE 9 tir')
surf(angle,position,curent1)
hold on 
surf(angle,position,curent2)
legend('I1 ANTIPHASE','I2 ANTIPHASE')
xlabel('angle')
ylabel('position')
%%
filename = 'data_10tir.xlsx';
start = 3;
power = zeros(11,10);
curent = zeros(11,10);
angle = 0:10:90;
position = -25:5:25;
for i=1:11
    rangex = 'D'+string(start+ ( i - 1 ) * 10 )+':D'+string(start + i * 10 );
    power(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'C'+string(start+ ( i - 1 ) * 10 )+':C'+string(start + i * 10 );
    curent(i,:) = table2array(readtable(filename,'Range',rangex));
end
figure('Name','single transmiter recieved power 10 tir')
surf(angle,position,power)
legend('Power')
xlabel('angle')
ylabel('position')
figure('Name','single transmiter current 10 tir')
surf(angle,position,curent)
legend('current')
xlabel('angle')
ylabel('position')