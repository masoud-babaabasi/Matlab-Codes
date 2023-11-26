clear 
filename = 'data_9tir.xlsx';
start = 3;
Vr_single = zeros(11,10);
It_single = zeros(11,10);
angle = 0:10:90;
position_single = -25:5:25;
for i=1:11
    rangex = 'D'+string(start+ ( i - 1 ) * 10 )+':D'+string(start + i * 10 );
    Vr_single(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'C'+string(start+ ( i - 1 ) * 10 )+':C'+string(start + i * 10 );
    It_single(i,:) = table2array(readtable(filename,'Range',rangex));
end
position_dual = 0:5:50;
Vr_IN = zeros(11,10);
I1_IN = zeros(11,10);
I2_IN = zeros(11,10);
for i=1:11
    rangex = 'L'+string(start+ ( i - 1 ) * 10 )+':L'+string(start + i * 10 );
    Vr_IN(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'J'+string(start+ ( i - 1 ) * 10 )+':J'+string(start + i * 10 );
    I1_IN(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'K'+string(start+ ( i - 1 ) * 10 )+':K'+string(start + i * 10 );
    I2_IN(i,:) = table2array(readtable(filename,'Range',rangex));
end

Vr_ANTI = zeros(11,10);
I1_ANTI = zeros(11,10);
I2_ANTI = zeros(11,10);
for i=1:11
    rangex = 'U'+string(start+ ( i - 1 ) * 10 )+':U'+string(start + i * 10 );
    Vr_ANTI(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'S'+string(start+ ( i - 1 ) * 10 )+':S'+string(start + i * 10 );
    I1_ANTI(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'T'+string(start+ ( i - 1 ) * 10 )+':T'+string(start + i * 10 );
    I2_ANTI(i,:) = table2array(readtable(filename,'Range',rangex));
end
%%
R = 25*sqrt(2);
w = 1;
space = 0.35;
Nt = 10;
%% transmiter coils
x0 = 0;
y0 = 0;
coil_t = zeros(2,3,Nt*4);
p=[[0,0];[53,0]];
for k= 1:2
    x0 = p(k,1) ;
    y0 = p(k,2) ;
    r = R;
    for i= 1:Nt
       for j= 1:4
          theta = (j-1)*2*pi/4+pi/4;
          coil_t(k,:,(i-1)*4+j) = [r*cos(theta)+x0; r*sin(theta)+y0;0];
       end
        r = r - 2*( w + space)*tan(pi/4);
    end
%     plot3(squeeze(coil_t(k,1,:)),squeeze(coil_t(k,2,:)),squeeze(coil_t(k,3,:)),'r')
%     hold on
end
% reciver coil 
Nr = 7;
R_rec = 15*sqrt(2);
xr = 25;
yr = 0;
zr = 40;
x_angle = 0;
y_angle = 45*pi/180;
z_angle = 0;
coil_r = zeros(3,Nr*4);
r = R_rec;
for i= 1:Nr
       for j= 1:4
          theta = (j-1)*2*pi/4+pi/4;
          coil_r(:,(i-1)*4+j) = [r*cos(theta)+xr; r*sin(theta)+yr;zr];
       end
        r = r - 2*( w + space)*tan(pi/4);
end
ds = [0;0;1]*R_rec;
r_z = [cos(z_angle) -sin(z_angle) 0 ; sin(z_angle) cos(z_angle) 0 ;0 0 1];
r_y = [cos(y_angle) 0 sin(y_angle)  ; 0 1 0 ; -sin(y_angle) 0 cos(y_angle)];
r_x = [1 0 0; 0 cos(x_angle) -sin(x_angle);0 sin(x_angle) cos(x_angle) ];
r = r_x*r_y*r_z;
coil_r = (r*(coil_r-[xr;yr;zr]))+[xr;yr;zr]; 
ds =  r * ds;
%%
s = serialport('COM5',115200);
ts = 50*1e-3;
Twindow = 5;
t = 0:ts:Twindow;
n = length(t);
data = zeros(n,3);
i = 1;
y_max = 500;
figure('Name','Pot values');
 for j=1:4
        subplot(4,1,j);
        if(j~=4)
            h(j) = plot(t,data(:,j));
        end
        if(j < 3 )
            txt(j)=text(0,125,['adc = ',num2str(data(n,j))]);
            xlim([0 Twindow]);
            xticks([0:1:Twindow])
            ylim([0 y_max]);
            ylabel("coil"+string(j),'fontsize',14);
            grid on;
            %set(get(gca,'ylabel'),'rotation',0)
            yticks([0:100:y_max]);
        elseif(j == 3)
           txt(j)=text(0,125,['Vr = ',num2str(data(n,j))]);
            xlim([0 Twindow]);
            xticks([0:1:Twindow])
            ylim([0 2500]);
            ylabel("V reciever",'fontsize',14);
            grid on;
            %set(get(gca,'ylabel'),'rotation',0)
            yticks([0:500:2500]); 
        else
            for k= 1:2
                plot3(squeeze(coil_t(k,1,:)),squeeze(coil_t(k,2,:)),squeeze(coil_t(k,3,:)),'r')
                hold on
            end
            reciever = plot3(coil_r(1,:),coil_r(2,:),coil_r(3,:),'b');
            ds_plot = quiver3(xr,yr,zr,ds(1),ds(2),ds(3),'G');
            axis equal
        end
 end
drawnow;
%%
state = 'I';
write(s,double(state),"uint8");
read_data = zeros(1,6);
while(1)
    read_data = read(s,6,"uint8");
    read_data_6 = read_data(1:2:5) + read_data(2:2:6)*256;
    if(i <= n)
        data(n-i+1,:) = read_data_6;
    else
        data(1:n-1,:) = data(2:n,:);
        data(n,:) = read_data_6;
    end
    for j=1:3
        h(j).YData = data(:,j);
        if( j ~= 3 )
            txt(j).String = ['mean adc = ',num2str(round(mean(data(:,j))))];
        else
            txt(j).String = ['mean adc = ',num2str(round(mean(data(:,j)))) , ' , Vr = ',num2str(round(mean(data(:,j)))/4096*3.3/0.1568) ,' x= ',num2str(xr) ,', angle= ',num2str(round(y_angle*180/pi)) ];
        end        
    end
    i = i + 1;
    if( mod(i,100) == 0 )
        I1 = round(mean(data(:,1)));
        I2 = round(mean(data(:,2)));
        Vr = round(mean(data(:,3)));
        e1 = abs(I1_IN - I1).^2 + abs(I2_IN - I2).^2 +(Vr - Vr_IN).^2;
        [Min,Index] = min(e1);
        [Min2,Index2] = min(Min);
        coil_r = coil_r-[xr;yr;zr];
        xr = position_dual(Index(Index2));
%         coil_r = coil_r+[xr;yr;zr];
        ds = [0;0;1]*R_rec;
        coil_r = ( r \ coil_r );
        y_angle = -angle(Index2)*pi/180;
        r_y = [cos(y_angle) 0 sin(y_angle)  ; 0 1 0 ; -sin(y_angle) 0 cos(y_angle)];
        r = r_x*r_y*r_z;
        coil_r = (r * coil_r ) +[xr;yr;zr];
        ds = r * ds;
        reciever.XData = coil_r(1,:);
        reciever.YData = coil_r(2,:);
        reciever.ZData = coil_r(3,:);

        ds_plot.UData = ds(1);
        ds_plot.VData = ds(2);
        ds_plot.WData = ds(3);
        ds_plot.XData = xr;
        ds_plot.YData = yr;
        ds_plot.ZData = zr;
    end
end

