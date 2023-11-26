clear 
load('matris.mat')
R = 25e-3 ; %Meter
w = 1e-3;
z = squeeze(max(max(max(K))))*100;
C1 = [ 0 ; (R+w)/cos(pi/6) ; z ];
C2 = [(R+w) ; -(R+w)*tan(pi/6) ; z ];
C3 = [-(R+w); -(R+w)*tan(pi/6) ; z ];
t = 0:1e-3:2*pi;
coil = [R*cos(t); R*sin(t);zeros(1,length(t))];
Range_tilt = 0:pi/6:pi/2;
N_tilt = length(Range_tilt);
tilt = [repmat(Range_tilt,1,N_tilt) ; kron(Range_tilt,zeros(1,N_tilt))];

%%
RangeD = -3*R:5e-3:3*R;
Nd = length(RangeD);

%%
alpha = 0.3;
for n_k = 1:4
    figure('Name','coil'+string(n_k))
    plot3(coil(1,:)+C1(1),coil(2,:)+C1(2),coil(3,:)+C1(3),'m')
    hold on
    plot3(coil(1,:)+C2(1),coil(2,:)+C2(2),coil(3,:)+C2(3),'m')
    plot3(coil(1,:)+C3(1),coil(2,:)+C3(2),coil(3,:)+C3(3),'m')
    
    x_angle = tilt(1,n_k);
    y_angle = tilt(2,n_k);
    z_angle = 0;
    r_z = [cos(z_angle) -sin(z_angle) 0 ; sin(z_angle) cos(z_angle) 0 ;0 0 1];
    r_y = [cos(y_angle) 0 sin(y_angle)  ; 0 1 0 ; -sin(y_angle) 0 cos(y_angle)];
    r_x = [1 0 0; 0 cos(x_angle) -sin(x_angle);0 sin(x_angle) cos(x_angle) ];
    r = r_x*r_y*r_z;    
    coilr = r * coil;
    coilr = [coilr(1,:); coilr(2,:);coilr(3,:)+ones(1,length(t))*2*C2(3)];
%     plot3(coilr(1,:),coilr(2,:),coilr(3,:),'r')

    CO(:,:,1) = zeros(Nd); % red
    CO(:,:,2) = zeros(Nd); % green
    CO(:,:,3) = ones(Nd); % blue
    surf(RangeD+C1(1),RangeD+C1(2),squeeze((K(n_k,:,:)))'*100,CO,'FaceColor','interp','FaceAlpha',alpha);
    CO(:,:,1) = zeros(Nd); % red
    CO(:,:,3) = zeros(Nd); % green
    CO(:,:,2) = ones(Nd); % blue
    surf(RangeD+C2(1),RangeD+C2(2),squeeze((K(n_k,:,:)))'*100,CO,'FaceColor','interp','FaceAlpha',alpha);
    CO(:,:,3) = zeros(Nd); % red
    CO(:,:,2) = zeros(Nd); % green
    CO(:,:,1) = ones(Nd); % blue
    surf(RangeD+C3(1),RangeD+C3(2),squeeze((K(n_k,:,:)))'*100,CO,'FaceColor','interp','FaceAlpha',alpha);
end