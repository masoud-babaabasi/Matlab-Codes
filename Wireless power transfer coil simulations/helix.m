clear 
R = 25e-3;
N = 10;
w = 1e-3;
s = 2e-3;
l = N*(w+s);
I = 1;
dx = 10e-3;
Range = -2*l:dx:2*l;
Nr = length(Range);
P = [repmat(Range,1,Nr*Nr) ; kron(repmat(Range,1,Nr),ones(1,Nr)) ; kron(Range,ones(1,Nr*Nr))];
% t=0:5e-3:2*pi*N;
% Helix = [R*cos(t) ; R*sin(t) ; -l/2+t/2/pi*(w+s)];
% plot3(Helix(1,:),Helix(2,:),Helix(3,:))
C=[0 0 0];
B = zeros(size(P));
%%
for i=1:Nr
   for j =1:Nr
      for k = 1:Nr
          a = (k-1) + (j-1)*Nr + (i-1)*Nr*Nr + 1;
          if abs(sqrt((P(1,a)^2 + P(2,a)^2))-R) > (2e-3)
            B(:,a) = integral(@(theta)intgral_func4(theta,R,C,w,s,l,P(:,a),I),0,2*pi*N,'ArrayValued',true);
          end
      end
   end
end
%%
figure('Name','Helix coil')
quiver3(P(1,:),P(2,:),P(3,:),B(1,:),B(2,:),B(3,:),'B')
hold on
t=0:5e-3:2*pi*N;
Helix = [R*cos(t) ; R*sin(t) ; -l/2+t/2/pi*(w+s)];
plot3(Helix(1,:),Helix(2,:),Helix(3,:),'r')
axis equal