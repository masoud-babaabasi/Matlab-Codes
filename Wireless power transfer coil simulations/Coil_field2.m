clear 
R = 25e-3 ; %Meter
I = 1;
dx = 10e-3;
Range = -2*R:dx:2*R;
N = length(Range);
P = [repmat(Range,1,N*N) ; kron(repmat(Range,1,N),ones(1,N)) ; kron(Range,ones(1,N*N))];
B1 = zeros(size(P));
B2 = zeros(size(P));
B3 = zeros(size(P));
C1 = [sqrt(3)*R - R /2 ; 0 ; 0 ];
C2 = [ - R /2 ; -R ; 0 ];
C3 = [ - R /2 ; R ; 0 ];
%%
r = 5e-3;
for i=1:N
   for j =1:N
      for k = 1:N
          a = (k-1) + (j-1)*N + (i-1)*N*N + 1;
%           if abs(P(3,a)) >= 15e-3  
          if ((P(1,a)-C1(1))^2 + (P(2,a)-C1(2))^2 + P(3,a)^2 + R^2 - 2*R*sqrt((P(1,a)-C1(1))^2 + (P(2,a)-C1(2))^2)) > (r)^2 && ((P(1,a)-C2(1))^2 + (P(2,a)-C2(2))^2 + P(3,a)^2 + R^2 - 2*R*sqrt((P(1,a)-C2(1))^2 + (P(2,a)-C2(2))^2)) > (r)^2 && ((P(1,a)-C3(1))^2 + (P(2,a)-C3(2))^2 + P(3,a)^2 + R^2 - 2*R*sqrt((P(1,a)-C3(1))^2 + (P(2,a)-C3(2))^2)) > (r)^2
            B1(:,a) = integral(@(theta)intgral_func(theta,R,C1,P(:,a),I),0,2*pi,'ArrayValued',true);
            B2(:,a) = integral(@(theta)intgral_func(theta,R,C2,P(:,a),I),0,2*pi,'ArrayValued',true);
            B3(:,a) = integral(@(theta)intgral_func(theta,R,C3,P(:,a),I),2*pi,0,'ArrayValued',true);
            B32(:,a) = integral(@(theta)intgral_func(theta,R,C3,P(:,a),I),0,2*pi,'ArrayValued',true);
          end
      end
   end
end
B = B1 + B2 + B3 ;
B_2 = B1 + B2 + B32 ;
%%
B_unit = B ./ vecnorm(B);
B_norm = abs(log(vecnorm(B)+1));
B_norm( isinf(B_norm)) = 0;
B_log = B_unit .* B_norm;

B2_unit = B_2 ./ vecnorm(B_2);
B2_norm = abs(log(vecnorm(B_2)+1));
B2_norm( isinf(B2_norm)) = 0;
B2_log = B2_unit .* B2_norm;
figure('Name','Three coil')
quiver3(P(1,:),P(2,:),P(3,:),B_log(1,:),B_log(2,:),B_log(3,:),'B')
hold on
quiver3(P(1,:),P(2,:),P(3,:),B2_log(1,:),B2_log(2,:),B2_log(3,:),'r')
t = 0:1e-3:2*pi;
coil1 = [R*cos(t)+C1(1); R*sin(t)+C1(2);ones(1,length(t))*(C1(3))];
coil2 = [R*cos(t)+C2(1); R*sin(t)+C2(2);ones(1,length(t))*(C2(3))];
coil3 = [R*cos(t)+C3(1); R*sin(t)+C3(2);ones(1,length(t))*(C3(3))];
plot3(coil1(1,:),coil1(2,:),coil1(3,:))
plot3(coil2(1,:),coil2(2,:),coil2(3,:))
plot3(coil3(1,:),coil3(2,:),coil3(3,:))
dl1 = I*[-sin(t(1));cos(t(1));0]*R;
dl2 = -dl1;
quiver3(coil1(1),coil1(2),coil1(3),dl1(1),dl1(2),dl1(3),'G')
quiver3(coil2(1),coil2(2),coil2(3),dl1(1),dl1(2),dl1(3),'G')
quiver3(coil3(1),coil3(2),coil3(3),dl2(1),dl2(2),dl2(3),'G')