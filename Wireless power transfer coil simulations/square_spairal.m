clear 
L = 50e-3;
R = L / 2;
C = [0,0,0];
N =13;
w = 1e-3;
s = 0.35e-3;
dx = 10e-3;
I =1;
Range = -2*R:dx:2*R;
Nr = length(Range);
P = [repmat(Range,1,Nr*Nr) ; kron(repmat(Range,1,Nr),ones(1,Nr)) ; kron(Range,ones(1,Nr*Nr))];
B1 = zeros(size(P));
B2 = zeros(size(P));
B3 = zeros(size(P));
B4 = zeros(size(P));
%%
for i=1:Nr
   for j =1:Nr
       for h=1:Nr
          a = (h-1) + (j-1)*Nr + (i-1)*Nr*Nr + 1;
          for k = 0:N-1
              X1 = [ C(1)-R+(w+s)*(k-1) ; C(1)+R-(w+s)*k  ; 0 ];
              X2 = [ C(1)+R-(w+s)*k ; C(1)+R-(w+s)*k  ; 0 ];
              X3 = [ C(1)+R-(w+s)*k ; C(1)-R+(w+s)*k  ; 0 ];
              X4 = [ C(1)-R+(w+s)*k ; C(1)-R+(w+s)*k  ; 0 ];
              X5 = [ C(1)-R+(w+s)*k ; C(1)+R-(w+s)*(k+1) ; 0 ];
              con1 =  X1(1) < (P(1,a)-C(1)) && P(1,a)-C(1) < X2(1) + (w/2) && P(2,a)-C(2) < X1(2) + (w/2) && P(2,a)-C(2) > X1(2) - (w/2);
              con2 =  X2(1) - (w/2) < (P(1,a)-C(1)) && P(1,a)-C(1) < X2(1) + (w/2) && P(2,a)-C(2) < X2(2) + (w/2) && P(2,a)-C(2) > X3(2) - (w/2);
              con3 =  X3(1) + (w/2) > (P(1,a)-C(1)) && P(1,a)-C(1) > X4(1) - (w/2) && P(2,a)-C(2) < X3(2) + (w/2) && P(2,a)-C(2) > X3(2) - (w/2);
              con4 =  X4(1) - (w/2) < (P(1,a)-C(1)) && P(1,a)-C(1) < X4(1) + (w/2) && P(2,a)-C(2) < X5(2) + (w/2) && P(2,a)-C(2) > X4(2) - (w/2);
              if((con1 || con2 || con3 || con4) && P(3,a) == 0)
                  B1(:,a) = 0;
                  B2(:,a) = 0;
                  B3(:,a) = 0;
                  B4(:,a) = 0;
                  break;
              else
                  B1(:,a) = B1(:,a) + integral(@(t)intgral_func2(t,X1,X2,P(:,a),I),0,1,'ArrayValued',true);
                  B2(:,a) = B2(:,a) + integral(@(t)intgral_func2(t,X2,X3,P(:,a),I),0,1,'ArrayValued',true);
                  B3(:,a) = B3(:,a) + integral(@(t)intgral_func2(t,X3,X4,P(:,a),I),0,1,'ArrayValued',true);
                  B4(:,a) = B4(:,a) + integral(@(t)intgral_func2(t,X4,X5,P(:,a),I),0,1,'ArrayValued',true);
              end
         end
       end
   end
end
B = B1 + B2 + B3 + B4;
%%
figure('Name','line Field')
quiver3(P(1,:),P(2,:),P(3,:),B(1,:),B(2,:),B(3,:),'B')
hold on
t1 = 0:1e-3:1;
for k = 0:N-1
X1 = [ C(1)-R+(w+s)*(k-1) ; C(1)+R-(w+s)*k  ; 0 ];
X2 = [ C(1)+R-(w+s)*k ; C(1)+R-(w+s)*k  ; 0 ];
X3 = [ C(1)+R-(w+s)*k ; C(1)-R+(w+s)*k  ; 0 ];
X4 = [ C(1)-R+(w+s)*k ; C(1)-R+(w+s)*k  ; 0 ];
X5 = [ C(1)-R+(w+s)*k ; C(1)+R-(w+s)*(k+1) ; 0 ]; 
dl = [X2(1)-X1(1);X2(2)-X1(2);X2(3)-X1(3)];
line = [dl(1)*t1+X1(1) ; dl(2)*t1+X1(2) ; dl(3)*t1+X1(3)];
plot3(line(1,:),line(2,:),line(3,:),'R');

dl = [X3(1)-X2(1);X3(2)-X2(2);X3(3)-X2(3)];
line = [dl(1)*t1+X2(1) ; dl(2)*t1+X2(2) ; dl(3)*t1+X2(3)];
plot3(line(1,:),line(2,:),line(3,:),'R');

dl = [X4(1)-X3(1);X4(2)-X3(2);X4(3)-X3(3)];
line = [dl(1)*t1+X3(1) ; dl(2)*t1+X3(2) ; dl(3)*t1+X3(3)];
plot3(line(1,:),line(2,:),line(3,:),'R');

dl = [X5(1)-X4(1);X5(2)-X4(2);X5(3)-X4(3)];
line = [dl(1)*t1+X4(1) ; dl(2)*t1+X4(2) ; dl(3)*t1+X4(3)];
plot3(line(1,:),line(2,:),line(3,:),'R');
end