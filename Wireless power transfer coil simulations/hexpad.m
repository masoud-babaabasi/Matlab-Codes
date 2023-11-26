clear
R = 25e-3;
w = 1e-3;
s = 0.35e-3;
Nt = 12;
%% transmiter coils
x0 = 0;
y0 = 0;
coil_t = zeros(7,3,Nt*6);
for k= 1:7
    r = R;
    for i= 1:Nt
       for j= 1:6
          theta = (j-1)*2*pi/6;
          coil_t(k,:,(i-1)*6+j) = [r*cos(theta)+x0; r*sin(theta)+y0;0];
       end
        r = r - 2*( w + s)*tan(pi/6);
    end
    plot3(squeeze(coil_t(k,1,:)),squeeze(coil_t(k,2,:)),squeeze(coil_t(k,3,:)),'r')
    hold on
    x0 = cos(k * 2 * pi / 6 + pi / 6) * (2 * R*cos(pi/6)+ 3e-3) ;
    y0 = sin(k * 2 * pi / 6 + pi / 6) * (2 * R*cos(pi/6)+ 3e-3) ;
end
%% reciver coil 
Nr_poly = 4;
Nr = 7;
R_rec = 15e-3*sqrt(2);
xr = 0;
yr = 0;
zr = 30e-3;
x_angle = 0;
y_angle = 0;
z_angle = 0;
coil_r = zeros(3,Nr*(Nr_poly+1));
r = R_rec;
for i=1:Nr
   for j=1:Nr_poly
       theta = (j-1)*2*pi/Nr_poly;
       coil_r(:,(i-1)*(Nr_poly+1)+j) = [r*cos(theta)+xr; r*sin(theta)+yr;zr];
   end
   coil_r(:,i*(Nr_poly+1)) = [r*cos(theta+2*pi/Nr_poly)+xr-(w+s)/sin(pi*(Nr_poly-2)/2/Nr_poly)/2; r*sin(theta+2*pi/Nr_poly)+yr-(w+s)/cos(pi*(Nr_poly-2)/2/Nr_poly)/2;zr];
   r = r - ( w + s)/(cos(pi/Nr_poly));
end
coil_r = [cos(pi/4) -sin(pi/4) 0 ; sin(pi/4) cos(pi/4) 0 ;0 0 1]*coil_r;
ds = [0;0;1]*R_rec;
reciever = plot3(coil_r(1,:),coil_r(2,:),coil_r(3,:),'b');
ds_plot = quiver3(xr,yr,zr,ds(1),ds(2),ds(3),'G');
axis equal
%% rotation
prompt = {'X(mm)','Y(mm)','Z(mm)','x_angel(d)','y_angel(d)','z_angel(d)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {num2str(xr*1e3),num2str(yr*1e3),num2str(zr*1e3),num2str(x_angle*180/pi),num2str(y_angle*180/pi),num2str(z_angle*180/pi)};
answer = inputdlg(prompt,dlgtitle,dims,definput);
xr2 = str2num(answer{1})*1e-3;
yr2 = str2num(answer{2})*1e-3;
zr2 = str2num(answer{3})*1e-3;
coil_r = coil_r + [xr2-xr ;yr2-yr;zr2-zr];
xr = xr2;
yr = yr2;
zr = zr2;
clear xr2 yr2 zr2

r_z = [cos(z_angle) -sin(z_angle) 0 ; sin(z_angle) cos(z_angle) 0 ;0 0 1];
r_y = [cos(y_angle) 0 sin(y_angle)  ; 0 1 0 ; -sin(y_angle) 0 cos(y_angle)];
r_x = [1 0 0; 0 cos(x_angle) -sin(x_angle);0 sin(x_angle) cos(x_angle) ];
r = r_x*r_y*r_z;
coil_r = (r\(coil_r-[xr;yr;zr]))+[xr;yr;zr]; 
ds =  r \ ds;

x_angle = str2num(answer{4})*pi/180;
y_angle = str2num(answer{5})*pi/180;
z_angle = str2num(answer{6})*pi/180;
r_z = [cos(z_angle) -sin(z_angle) 0 ; sin(z_angle) cos(z_angle) 0 ;0 0 1];
r_y = [cos(y_angle) 0 sin(y_angle)  ; 0 1 0 ; -sin(y_angle) 0 cos(y_angle)];
r_x = [1 0 0; 0 cos(x_angle) -sin(x_angle);0 sin(x_angle) cos(x_angle) ];
coil_r = (r_x * r_y * r_z*(coil_r-[xr;yr;zr]))+[xr;yr;zr]; 

ds = r_x*r_y*r_z * ds;
reciever.XData = coil_r(1,:);
reciever.YData = coil_r(2,:);
reciever.ZData = coil_r(3,:);

ds_plot.UData = ds(1);
ds_plot.VData = ds(2);
ds_plot.WData = ds(3);
ds_plot.XData = xr;
ds_plot.YData = yr;
ds_plot.ZData = zr;