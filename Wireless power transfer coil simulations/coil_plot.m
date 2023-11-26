clear 
figure
R1 = 50e-3/2 ; %Meter
w = 1e-3;
s = 0.35e-3;
C1 = [ 0 ; 0 ; 0 ];
R2 = 30e-3/2 ; %Meter
C2 = [ 0; 0 ; 40e-3 ];
n1 = 10;
n2 = 5;
t = 0:1e-3:2*pi*n1;
R = R1 - (n1 * ( w + s ) ) + t / 2 / pi * ( w + s );
dl = [0 , 1 , 0];
coil = [R.*cos(t)+C1(1); R.*sin(t)+C1(2);ones(1,length(t))*(C1(3))];
    coil = coil / R1;
    plot3(coil(1,:),coil(2,:),coil(3,:),'b')
    hold on 
    plot3(coil(1,:),coil(2,:)+2,coil(3,:),'b')
    x_angle = pi/2;
    y_angle = 0;
    z_angle = 0;
    r_z = [cos(z_angle) -sin(z_angle) 0 ; sin(z_angle) cos(z_angle) 0 ;0 0 1];
    r_y = [cos(y_angle) 0 sin(y_angle)  ; 0 1 0 ; -sin(y_angle) 0 cos(y_angle)];
    r_x = [1 0 0; 0 cos(x_angle) -sin(x_angle);0 sin(x_angle) cos(x_angle) ];
    r = r_x*r_y*r_z;  
    t2 = 0:1e-3:2*pi*n2;
    R = R2 - (n2 * ( w + s ) ) + t2 / 2 / pi * ( w + s );
    coil = [R.*cos(t2)+C1(1); R.*sin(t2)+C1(2);ones(1,length(t2))*(C1(3))];
    coil = coil / R1;
    coil = r * coil;
    coil = [coil(1,:); coil(2,:);coil(3,:)+ones(1,length(t2))*C2(3)/2/R1];
    plot3(coil(1,:),coil(2,:),coil(3,:),'r')
    axis equal
    xticks([0])
    yticks([0 2]) 
    yticklabels(["0" , "2R"])
    zticks([0 C2(3)/2/R1])
    zticklabels(["0" , string(C2(3)/R1)+"R"])
    xlabel('X')
    ylabel('Y')
%     quiver3(C2(1)/2/R1,C2(2)/2/R1,C2(3)/2/R1,dl(1),dl(2),dl(3),'G')
    grid on 
%%
R1 = 50e-3/2 ; %Meter
w = 1e-3;
s = 0.35e-3;
C1 = [ 0 ; 0 ; 0 ];
R2 = 30e-3/2 ; %Meter
C2 = [ 0; 0 ; 40e-3 ];
n1 = 10;
n2 = 5;
t = 0:1e-3:2*pi*n1;
R = R1 - (n1 * ( w + s ) ) + t / 2 / pi * ( w + s );
dl = [0 , 1 , 0];
coil = [R.*cos(t)+C1(1); R.*sin(t)+C1(2);ones(1,length(t))*(C1(3))];
coil = coil / R1;
plot(coil(1,:),coil(2,:),'b')
hold on 
plot(coil(1,:)+1.6,coil(2,:),'r')
axis equal
xticks([0 1.5]) 
xticklabels(["0" , "1.5R"])

yticks([-1:1]) 
yticklabels(["1R" , "0" , "1R"])
grid on