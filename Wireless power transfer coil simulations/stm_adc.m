clear 
prompt = {'Serial port:','baudrate:','sampling time(ms)','time window(s)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'COM5','115200','50','5'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
s = serialport(string(answer{1}),str2num(answer{2}));
ts = str2num(answer{3})*1e-3;
Twindow = str2num(answer{4});
t = 0:ts:Twindow;
n = length(t);
data = zeros(n,8);
i = 1;
y_max = 2200;
figure('Name','Pot values');
 for j=1:8
        scrollsubplot(8,1,[3*j-2,3*j]);
        h(j) = plot(t,data(:,j));
        if(j ~= 8 )
            txt(j)=text(0,125,['adc = ',num2str(data(n,j))]);
            xlim([0 Twindow]);
            xticks([0:1:Twindow])
            ylim([0 y_max]);
            ylabel("coil"+string(j),'fontsize',14);
            grid on;
            %set(get(gca,'ylabel'),'rotation',0)
            yticks([0:200:y_max]);
        else
           txt(j)=text(0,125,['Vr = ',num2str(data(n,j))]);
            xlim([0 Twindow]);
            xticks([0:1:Twindow])
            ylim([0 2000]);
            ylabel("V reciever",'fontsize',14);
            grid on;
            %set(get(gca,'ylabel'),'rotation',0)
            yticks([0:200:2000]); 
        end
 end
drawnow;
%%
read_data = zeros(1,14);
read_data_16 = zeros(1,7);
while(1)
    read_data = read(s,14,"uint8");
    read_data_16 = read_data(1:2:13) + read_data(2:2:14)*256;
    data(1:n-1,:) = data(2:n,:);
    data(n,1:7) = read_data_16;

    for j=1:8
        if( j ~= 8 )
            h(j).YData = data(:,j);
            txt(j).String = ['mean adc = ',num2str(round(mean(data(:,j))))];
        else
            txt(j).String = ['mean adc = ',num2str(round(mean(data(:,j)))) , ' , Vr = ',num2str(round(mean(data(:,j)))/4096*3.3/0.1568)];
        end        
    end
    i = i + 1;
end

