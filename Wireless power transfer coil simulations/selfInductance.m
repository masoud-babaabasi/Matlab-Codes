clear 
R = 25e-3 ; %Meter
I = 1;
dx = 0.5e-3;
w = 1e-3;
Range = -(R-w/2):dx:(R-w/2);
N = length(Range);
P = [repmat(Range,1,N) ; kron(Range,ones(1,N)) ; zeros(1,N*N)];
B1 = zeros(size(P));
B2 = zeros(size(P));
B3 = zeros(size(P));
B4 = zeros(size(P));
X1 = [ -R ; -R ; 0 ];
X2 = [ -R ;  R ; 0 ];
X3 = [ R ;  R ; 0 ];
X4 = [ R ; -R ; 0 ];
%%
for i=1:N
   for j =1:N
          a = (j-1) + (i-1)*N + 1;
          B1(:,a) = integral(@(t)intgral_func2(t,X1,X2,P(:,a),I),0,1,'ArrayValued',true);
          B2(:,a) = integral(@(t)intgral_func2(t,X2,X3,P(:,a),I),0,1,'ArrayValued',true);
          B3(:,a) = integral(@(t)intgral_func2(t,X3,X4,P(:,a),I),0,1,'ArrayValued',true);
          B4(:,a) = integral(@(t)intgral_func2(t,X4,X1,P(:,a),I),0,1,'ArrayValued',true);
   end
end
%%
B = B1 + B2 + B3 + B4;
L = abs(sum(B(3,:))) * dx * dx / I;
figure('Name','line Field')
quiver3(P(1,:),P(2,:),P(3,:),B(1,:),B(2,:),B(3,:),'B')
hold on
t1 = 0:1e-3:1;
dl = [X2(1)-X1(1);X2(2)-X1(2);X2(3)-X1(3)];
line = [dl(1)*t1+X1(1) ; dl(2)*t1+X1(2) ; dl(3)*t1+X1(3)];
plot3(line(1,:),line(2,:),line(3,:),'R');

dl = [X3(1)-X2(1);X3(2)-X2(2);X3(3)-X2(3)];
line = [dl(1)*t1+X2(1) ; dl(2)*t1+X2(2) ; dl(3)*t1+X2(3)];
plot3(line(1,:),line(2,:),line(3,:),'R');

dl = [X4(1)-X3(1);X4(2)-X3(2);X4(3)-X3(3)];
line = [dl(1)*t1+X3(1) ; dl(2)*t1+X3(2) ; dl(3)*t1+X3(3)];
plot3(line(1,:),line(2,:),line(3,:),'R');

dl = [X1(1)-X4(1);X1(2)-X1(2);X1(3)-X4(3)];
line = [dl(1)*t1+X4(1) ; dl(2)*t1+X4(2) ; dl(3)*t1+X4(3)];
plot3(line(1,:),line(2,:),line(3,:),'R');