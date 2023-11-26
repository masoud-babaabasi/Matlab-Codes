clear 
clear 
% list = blelist
bluetooth_Name = "reciever coil";
b = ble(bluetooth_Name);
charac = b.Characteristics;
ServUUID = charac(5,2);
CharUUID = charac(5,4);
c = characteristic(b,ServUUID{:,:},CharUUID{:,:});
c2 = characteristic(b,"1801","2A05");
% subscribe(c2)

ts = 50*1e-3;
Twindow = 5;
t = 0:ts:Twindow;
n = length(t);
data = zeros(n,3);
i = 1;
y_max = 500;
figure('Name','Pot values');
 for j=1:3
        subplot(3,1,j);
        h(j) = plot(t,data(:,j));
        if(j == 1 )
            txt(j)=text(0,5,['Vin = ',num2str(data(n,j))]);
            xlim([0 Twindow]);
            xticks([0:1:Twindow])
            ylim([0 16]);
            ylabel("Vin",'fontsize',14);
            grid on;
            yticks([0:2:16]);
        elseif(j == 2)
           txt(j)=text(0,20,['Iin = ',num2str(data(n,j))]);
            xlim([0 Twindow]);
            xticks([0:1:Twindow])
            ylim([0 50]);
            ylabel("Iin",'fontsize',14);
            grid on;
            yticks([0:5:50]); 
         elseif(j == 3)
           txt(j)=text(0,125,['Power = ',num2str(data(n,j))]);
            xlim([0 Twindow]);
            xticks([0:1:Twindow])
            ylim([0 350]);
            ylabel("Power",'fontsize',14);
            grid on;
            yticks([0:50:300]); 
        end
 end
drawnow;
data_in = zeros(1,2);
while(1)
    ble_data = read(c);
%     if length(ble_data ) ~= 2 
%         continue;
%     end
%     char(data)
%     data_in =sscanf(string(char(ble_data)),"vin=%f , Iin = %fmA");
    data_in(1) = ble_data(1) * 14 / 255 + 3 ;
    data_in(2) = ble_data(2) * 60 / 255 + 10 ;
%     vin = data_in(1);
%     Iin = data_in(2);
    data(1:n-1,:) = data(2:n,:);
    data(n,1:2) = data_in;
    for j=1:2
        h(j).YData = data(:,j);
        if( j == 1 )
            txt(j).String = ['Vin  = ', num2str(data(n,j))];
        else
            txt(j).String = ['Iin = ', num2str(data(n,j)) ];
        end        
    end
        h(3).YData = data(:,1) .* data(:,2);
        txt(3).String = ['Power = ', num2str(data(n,1)*data(n,2)) ];
     i = i + 1;
     drawnow;
end
%write(c,data)