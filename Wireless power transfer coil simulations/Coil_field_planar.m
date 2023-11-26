clear 
R = 25e-3 ; %Meter
w = 1e-3;
s = 0.35e-3;
I = 1;
% dx = 10e-3;
dx = R/3;
err = 10e-3;

w = w / R;
s = s / R;
dx = dx / R;
err = err / R;
R = 1;
Range = -3*R:dx:3*R;
N = length(Range);
P = [zeros(1,N*N) ; repmat(Range,1,N) ; kron(Range,ones(1,N))];

%%
w = w / R;
s = s / R;
dx = dx / R;
R = 1;
B1_in = zeros(size(P));
B2_in = zeros(size(P));
B2_anti = zeros(size(P));
C1 = [ 0 ; -R ; 0 ];
C2 = [ 0 ; R ; 0 ];
for i=1:N
   for j =1:N
          a =  (j-1) + (i-1)*N + 1;
%           if ((abs(sqrt(P(1,a)^2 + P(2,a)^2) - R  ) >= 10e-3) || (abs(P(3,a)- C(3)) >= 10e-3 )) && (P(3,a) ~= 0 )
          if (P(1,a)^2 + (P(2,a)-C1(2))^2 + P(3,a)^2 + R^2 - 2*R*sqrt(P(1,a)^2 + (P(2,a)-C1(2))^2)) > (err)^2 && (P(1,a)^2 + (P(2,a)-C2(2))^2 + P(3,a)^2 + R^2 - 2*R*sqrt(P(1,a)^2 + (P(2,a)-C2(2))^2)) > (err)^2 
              B1_in(:,a) = integral(@(theta)intgral_func(theta,R,C1,P(:,a),I,eye(3)),0,2*pi,'ArrayValued',true);
              B2_in(:,a) = integral(@(theta)intgral_func(theta,R,C2,P(:,a),I,eye(3)),0,2*pi,'ArrayValued',true);
              B2_anti(:,a) = integral(@(theta)intgral_func(theta,R,C2,P(:,a),I,eye(3)),2*pi,0,'ArrayValued',true);
%           B_in(:,a) = log10(abs(B_in(:,a))) ./ sign(B_in(:,a));
%           B_in2(:,a) = log10(abs(B_in2(:,a))) ./ sign(B_in2(:,a));
%           B_anti(:,a) = log10(abs(B_anti(:,a))) ./ sign(B_anti(:,a));
          end
   end
end
%%
B_in = B1_in + B2_in;
B_anti = B1_in + B2_anti;
figure('Name','Two coil Antiphase')
quiver(P(2,:),P(3,:),B_anti(2,:),B_anti(3,:),'B')
hold on 
marker_size = 10;
plot(2*R,0,'o','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(2*R,0,'.','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(w+s,0,'o','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(w+s,0,'+','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)

plot(-(w+s),0,'o','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
plot(-(w+s),0,'+','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
plot(-2*R,0,'o','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
plot(-2*R,0,'.','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
xlim([-3 3])
ylim([-2 2])
xlabel('Y')
ylabel('Z')
figure('Name','Two coil Antiphase')
quiver(P(2,:),P(3,:),B_in(2,:),B_in(3,:),'B')
hold on 
plot(2*R,0,'o','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(2*R,0,'+','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(w+s,0,'o','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(w+s,0,'.','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)

plot(-(w+s),0,'o','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
plot(-(w+s),0,'+','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
plot(-2*R,0,'o','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
plot(-2*R,0,'.','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
xlim([-3 3])
ylim([-2 2])
xlabel('Y')
ylabel('Z')
%%
Range = -5*R:dx:5*R;
N = length(Range);
clear p
P = [zeros(1,N*N) ; repmat(Range,1,N) ; kron(Range,ones(1,N))];
B1_in = zeros(size(P));
B2_anti = zeros(size(P));
C1 = [ 0 ; -2*R ; 0 ];
C2 = [ 0 ; 2*R ; 0 ];
for i=1:N
   for j =1:N
          a =  (j-1) + (i-1)*N + 1;
%           if ((abs(sqrt(P(1,a)^2 + P(2,a)^2) - R  ) >= 10e-3) || (abs(P(3,a)- C(3)) >= 10e-3 )) && (P(3,a) ~= 0 )
          if (P(1,a)^2 + (P(2,a)-C1(2))^2 + P(3,a)^2 + R^2 - 2*R*sqrt(P(1,a)^2 + (P(2,a)-C1(2))^2)) > (err)^2 && (P(1,a)^2 + (P(2,a)-C2(2))^2 + P(3,a)^2 + R^2 - 2*R*sqrt(P(1,a)^2 + (P(2,a)-C2(2))^2)) > (err)^2 
              B1_in(:,a) = integral(@(theta)intgral_func(theta,R,C1,P(:,a),I,eye(3)),0,2*pi,'ArrayValued',true);
%               B2_in(:,a) = integral(@(theta)intgral_func(theta,R,C2,P(:,a),I,eye(3)),0,2*pi,'ArrayValued',true);
              B2_anti(:,a) = integral(@(theta)intgral_func(theta,R,C2,P(:,a),I,eye(3)),2*pi,0,'ArrayValued',true);
%           B_in(:,a) = log10(abs(B_in(:,a))) ./ sign(B_in(:,a));
%           B_in2(:,a) = log10(abs(B_in2(:,a))) ./ sign(B_in2(:,a));
%           B_anti(:,a) = log10(abs(B_anti(:,a))) ./ sign(B_anti(:,a));
          end
   end
end
%%
B_anti = B1_in + B2_anti;
figure('Name','Two coil Antiphase')
quiver(P(2,:),P(3,:),B_anti(2,:),B_anti(3,:),'B')
hold on 
plot(3*R,0,'o','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(3*R,0,'.','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(R+w+s,0,'o','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(R+w+s,0,'+','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)

plot(-R-w-s,0,'o','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
plot(-R-w-s,0,'+','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
plot(-3*R,0,'o','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
plot(-3*R,0,'.','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)

plot(-R+w+s,0,'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',marker_size)
plot(R-w-s,0,'o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',marker_size)
xlim([-4 4])
ylim([-2 2])
xlabel('Y')
ylabel('Z')