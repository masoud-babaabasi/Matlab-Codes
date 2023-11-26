clear 
R = 25e-3 ; %Meter
I = 1;
dx = 5e-3;
dy = 5e-3;
dz = 10e-3;
Range_x = -R:dx:R;
Range_y = -R:dy:R;
Range_z = -2*R:dz:2*R;
Nx = length(Range_x);
Ny = length(Range_y);
Nz = length(Range_z);
P = [repmat(Range_x,1,Ny*Nz) ; kron(repmat(Range_y,1,Nz),ones(1,Nx)) ; kron(Range_z,ones(1,Nx*Ny))];
B1 = zeros(size(P));
X1 = [ 0 ; 0 ; -R ];
X2 = [0 ; 0 ; R ];
%%
for i=1:Nz
   for j =1:Ny
      for k = 1:Nx
          a = (k-1) + (j-1)*Nx + (i-1)*Nx*Ny + 1;
          
          B1(:,a) = integral(@(t)intgral_func2(t,X1,X2,P(:,a),I),0,1,'ArrayValued',true);
      end
   end
end
%%
figure('Name','line Field')
quiver3(P(1,:),P(2,:),P(3,:),B1(1,:),B1(2,:),B1(3,:),'B')
hold on
t1 = 0:1e-3:1;
dl = [X2(1)-X1(1);X2(2)-X1(2);X2(3)-X1(3)];
% dl = dl / norm(dl);
line = [dl(1)*t1+X1(1) ; dl(2)*t1+X1(2) ; dl(3)*t1+X1(3)];
plot3(line(1,:),line(2,:),line(3,:));