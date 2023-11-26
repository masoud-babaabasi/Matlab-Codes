clear 
R = 25e-3 ; %Meter
w = 1e-3;
s = 0.35e-3;
I = 1;
% dx = 10e-3;
dx = R/3;
err = 10e-3;
H = 30e-3;

w = w / R;
s = s / R;
dx = dx / R;
err = err / R;
H = H / R;
R = 1;
Range = -4*R:dx:4*R;
N = length(Range);
P = [repmat(Range,1,N) ; kron(Range,ones(1,N)) ; H*ones(1,N*N)];

%%
w = w / R;
s = s / R;
dx = dx / R;
R = 1;
B1_in = zeros(size(P));
B2_in = zeros(size(P));
B3_in = zeros(size(P));
B3_anti = zeros(size(P));
C1 = [ - R /2 ; +R ; 0 ];
C2 = [ - R /2 ; -R ; 0 ];
C3 = [sqrt(3)*R - R /2 ; 0 ; 0 ];
for i=1:N
   for j =1:N
          a =  (j-1) + (i-1)*N + 1;
%           if ((abs(sqrt(P(1,a)^2 + P(2,a)^2) - R  ) >= 10e-3) || (abs(P(3,a)- C(3)) >= 10e-3 )) && (P(3,a) ~= 0 )
          if (P(1,a)^2 + (P(2,a)-C1(2))^2 + P(3,a)^2 + R^2 - 2*R*sqrt(P(1,a)^2 + (P(2,a)-C1(2))^2)) > (err)^2 && (P(1,a)^2 + (P(2,a)-C2(2))^2 + P(3,a)^2 + R^2 - 2*R*sqrt(P(1,a)^2 + (P(2,a)-C2(2))^2)) > (err)^2 
              B1_in(:,a) = integral(@(theta)intgral_func(theta,R,C1,P(:,a),I,eye(3)),0,2*pi,'ArrayValued',true);
              B2_in(:,a) = integral(@(theta)intgral_func(theta,R,C2,P(:,a),I,eye(3)),0,2*pi,'ArrayValued',true);
              B3_in(:,a) = integral(@(theta)intgral_func(theta,R,C3,P(:,a),I,eye(3)),0,2*pi,'ArrayValued',true);
              B3_anti(:,a) = integral(@(theta)intgral_func(theta,R,C3,P(:,a),I,eye(3)),2*pi,0,'ArrayValued',true);
%           B_in(:,a) = log10(abs(B_in(:,a))) ./ sign(B_in(:,a));
%           B_in2(:,a) = log10(abs(B_in2(:,a))) ./ sign(B_in2(:,a));
%           B_anti(:,a) = log10(abs(B_anti(:,a))) ./ sign(B_anti(:,a));
          end
   end
end
%%
B_in = B1_in + B2_in + B3_in;
B_anti = B1_in + B2_in + B3_anti;
figure('Name','Three coil Antiphase')
quiver(P(1,:),P(2,:),B_anti(1,:),B_anti(2,:),'B')
hold on 
t = 0:1e-3:2*pi;
coil1 = [R*cos(t)+C1(1); R*sin(t)+C1(2)];
coil2 = [R*cos(t)+C2(1); R*sin(t)+C2(2)];
coil3 = [R*cos(t)+C3(1); R*sin(t)+C3(2)];
plot(coil1(1,:),coil1(2,:),'B')
plot(coil2(1,:),coil2(2,:),'B')
plot(coil3(1,:),coil3(2,:),'R')
axis equal
xlim([-2 3])
ylim([-2.5 2.5])
marker_size = 8;
plot(0.5*R,1,'^','MarkerFaceColor','blue','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(0.5*R,-1,'^','MarkerFaceColor','blue','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(sqrt(3)*R + R /2 ,0,'v','MarkerFaceColor','red','MarkerEdgeColor','red','MarkerSize',marker_size)
xlabel('X')
ylabel('Y')
figure('Name','Three coil Inphase')
quiver(P(1,:),P(2,:),B_in(1,:),B_in(2,:),'B')
hold on 
t = 0:1e-3:2*pi;
coil1 = [R*cos(t)+C1(1); R*sin(t)+C1(2)];
coil2 = [R*cos(t)+C2(1); R*sin(t)+C2(2)];
coil3 = [R*cos(t)+C3(1); R*sin(t)+C3(2)];
plot(coil1(1,:),coil1(2,:),'B')
plot(coil2(1,:),coil2(2,:),'B')
plot(coil3(1,:),coil3(2,:),'B')
axis equal
xlim([-2 3])
ylim([-2.5 2.5])
marker_size = 8;
plot(0.5*R,1,'^','MarkerFaceColor','blue','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(0.5*R,-1,'^','MarkerFaceColor','blue','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(sqrt(3)*R + R /2 ,0,'^','MarkerFaceColor','blue','MarkerEdgeColor','blue','MarkerSize',marker_size)