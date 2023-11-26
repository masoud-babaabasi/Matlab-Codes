clear 
R = 25e-3 ; %Meter
I = 1;
dx = 10e-3;
Range = -2*R:dx:2*R;
N = length(Range);
P = [repmat(Range,1,N*N) ; kron(repmat(Range,1,N),ones(1,N)) ; kron(Range,ones(1,N*N))];
B1 = zeros(size(P));
C = [ 0 ; 0 ; 0 ]*1e-3;
%**********************
%%
x_angle = pi/2;
y_angle = pi/2;
z_angle = 0;
r_z = [cos(z_angle) -sin(z_angle) 0 ; sin(z_angle) cos(z_angle) 0 ;0 0 1];
r_y = [cos(y_angle) 0 sin(y_angle)  ; 0 1 0 ; -sin(y_angle) 0 cos(y_angle)];
r_x = [1 0 0; 0 cos(x_angle) -sin(x_angle);0 sin(x_angle) cos(x_angle) ];
r = r_x*r_y*r_z;
for i=1:N
   for j =1:N
      for k = 1:N
          a = (k-1) + (j-1)*N + (i-1)*N*N + 1;
%           if ((abs(sqrt(P(1,a)^2 + P(2,a)^2) - R  ) >= 15e-3) || (abs(P(3,a)- C(3)) >= 15e-3 )) && P(3,a) ~= 0 
          if (P(1,a)^2 + P(2,a)^2 + P(3,a)^2 + R^2 - 2*R*sqrt(P(1,a)^2 + P(2,a)^2)) > (5e-3)^2
            B1(:,a) = integral(@(theta)intgral_func(theta,R,C,P(:,a),I,r),0,2*pi,'ArrayValued',true);
          end
      end
   end
end
%%
figure('Name','One coil')
quiver3(P(1,:),P(2,:),P(3,:),B1(1,:),B1(2,:),B1(3,:),'B')
hold on
t = 0:1e-3:2*pi;
coil = [R*cos(t)+C(1); R*sin(t)+C(2);ones(1,length(t))*(C(3))];
coil = r * coil;
plot3(coil(1,:),coil(2,:),coil(3,:),'r')
dl = I*[-sin(t(1));cos(t(1));0]*R;
dl = r * dl;
quiver3(coil(1),coil(2),coil(3),dl(1),dl(2),dl(3),'G')
%%
B_in = zeros(size(P));
B_in2 = zeros(size(P));
B_anti = zeros(size(P));
C1 = [ 0 ; -3*R/2 ; 0 ];
C2 = [ 0 ; 3*R/2 ; 0 ];
for i=1:N
   for j =1:N
      for k = 1:N
          a = (k-1) + (j-1)*N + (i-1)*N*N + 1;
%           if ((abs(sqrt(P(1,a)^2 + P(2,a)^2) - R  ) >= 10e-3) || (abs(P(3,a)- C(3)) >= 10e-3 )) && (P(3,a) ~= 0 )
          if (P(1,a)^2 + (P(2,a)-C1(2))^2 + P(3,a)^2 + R^2 - 2*R*sqrt(P(1,a)^2 + (P(2,a)-C1(2))^2)) > (15e-3)^2 && (P(1,a)^2 + (P(2,a)-C2(2))^2 + P(3,a)^2 + R^2 - 2*R*sqrt(P(1,a)^2 + (P(2,a)-C2(2))^2)) > (15e-3)^2 
              B_in(:,a) = integral(@(theta)intgral_func(theta,R,C1,P(:,a),I,eye(3)),0,2*pi,'ArrayValued',true);
              B_in2(:,a) = integral(@(theta)intgral_func(theta,R,C2,P(:,a),I,eye(3)),0,2*pi,'ArrayValued',true);
              B_anti(:,a) = integral(@(theta)intgral_func(theta,R,C2,P(:,a),I,eye(3)),2*pi,0,'ArrayValued',true);
%           B_in(:,a) = log10(abs(B_in(:,a))) ./ sign(B_in(:,a));
%           B_in2(:,a) = log10(abs(B_in2(:,a))) ./ sign(B_in2(:,a));
%           B_anti(:,a) = log10(abs(B_anti(:,a))) ./ sign(B_anti(:,a));
          end
      end
   end
end
%%
B2_anti = B_in + B_anti;
B2a_unit = B2_anti ./ vecnorm(B2_anti);
B2a_norm = abs(log(vecnorm(B2_anti)+1));
B2a_norm( isinf(B2a_norm)) = 0;
B2a_log = B2a_unit .* B2a_norm;
figure('Name','Two coil Antiphase')
quiver3(P(1,:),P(2,:),P(3,:),B2a_log(1,:),B2a_log(2,:),B2a_log(3,:),'B')
hold on
t = 0:1e-3:2*pi;
coil1 = [R*cos(t)+C1(1); R*sin(t)+C1(2);ones(1,length(t))*(C1(3))];
coil2 = [R*cos(t)+C2(1); R*sin(t)+C2(2);ones(1,length(t))*(C2(3))];
plot3(coil1(1,:),coil1(2,:),coil1(3,:))
plot3(coil2(1,:),coil2(2,:),coil2(3,:))

%%
B2_in = B_in + B_in2;
B2i_unit = B2_in ./ vecnorm(B2_in);
B2i_norm = abs(log(vecnorm(B2_in)+1));
B2i_norm( isinf(B2i_norm)) = 0;
B2i_log = B2i_unit .* B2i_norm;
figure('Name','Two coil Inphase')
quiver3(P(1,:),P(2,:),P(3,:),B2i_log(1,:),B2i_log(2,:),B2i_log(3,:),'B')
hold on
plot3(coil1(1,:),coil1(2,:),coil1(3,:))
plot3(coil2(1,:),coil2(2,:),coil2(3,:))
%%
figure('Name','Two coil field')
quiver3(P(1,:),P(2,:),P(3,:),B2a_log(1,:),B2a_log(2,:),B2a_log(3,:),'B')
hold on
quiver3(P(1,:),P(2,:),P(3,:),B2i_log(1,:),B2i_log(2,:),B2i_log(3,:),'R')
plot3(coil1(1,:),coil1(2,:),coil1(3,:))
plot3(coil2(1,:),coil2(2,:),coil2(3,:))