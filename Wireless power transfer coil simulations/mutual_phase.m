clear 
R1 = 50e-3/2 ; %Meter
w = 1e-3;
I = 1;
C1 = [ 0 ; 0 ; 0 ];
dx1 = 0.05 * R1;% dS = 0.8%
Range1 = -(R1-w/2):dx1:(R1-w/2);
N1 = length(Range1);
P1 = [repmat(Range1,1,N1) ; kron(Range1,ones(1,N1)) ; C1(3)*ones(1,N1*N1)];
B1 = zeros(size(P1));

R2 = 30e-3/2 ; %Meter
C2 = [ 0; 0 ; 40e-3 ];
dx2 = 0.05 * R2;% dS = 0.8%
Range2 = -(R2-w/2):dx2:(R2-w/2);
N2 = length(Range2);
P2 = [repmat(Range2,1,N2) ; kron(Range2,ones(1,N2)) ; C2(3)*ones(1,N2*N2)];
B2 = zeros(size(P2));
%%
tic;
for i=1:N1
   for j =1:N1
      a = (j-1) + (i-1)*N1 + 1;
      if( sqrt((P1(1,a) - C1(1))^2 + (P1(2,a) - C1(2))^2) < (R1 - w/2) )
        B1(:,a) = integral(@(theta)intgral_func(theta,R1,C1,P1(:,a),I,eye(3)),0,2*pi,'ArrayValued',true);
      end
   end
end
L1 = abs(sum(B1(3,:))) * dx1 * dx1 / I;

for i=1:N2
   for j =1:N2
      a = (j-1) + (i-1)*N2 + 1;
      if( sqrt((P2(1,a) - C2(1))^2 + (P2(2,a) - C2(2))^2) < (R2 - w/2) )
        B2(:,a) = integral(@(theta)intgral_func(theta,R2,C2,P2(:,a),I,eye(3)),0,2*pi,'ArrayValued',true);
      end
   end
end
L2 = abs(sum(B2(3,:))) * dx2 * dx2 / I;
%%
RangeD = 0:R1/4:2*R1;
Nd = length(RangeD);
% D = [repmat(RangeD,1,Nd) ; kron(RangeD,ones(1,Nd)) ; C2(3)*ones(1,Nd*Nd)];
Range_tilt = 0:pi/6:pi/2;
N_tilt = length(Range_tilt);
% tilt = [repmat(Range_tilt,1,N_tilt) ; kron(Range_tilt,ones(1,N_tilt))];
M21 = zeros(N_tilt,Nd);
for ii = 1:N_tilt
    x_angle = 0;
    y_angle = Range_tilt(ii);
    z_angle = 0;
    r_z = [cos(z_angle) -sin(z_angle) 0 ; sin(z_angle) cos(z_angle) 0 ;0 0 1];
    r_y = [cos(y_angle) 0 sin(y_angle)  ; 0 1 0 ; -sin(y_angle) 0 cos(y_angle)];
    r_x = [1 0 0; 0 cos(x_angle) -sin(x_angle);0 sin(x_angle) cos(x_angle) ];
    r = r_x*r_y*r_z;
    for di = 1:Nd
           C2 = [ RangeD(di); 0 ; 40e-3];
           B21 = zeros(size(P1));
           for i=1:N1
               for j =1:N1
                  a = (j-1) + (i-1)*N1 + 1;
                  P_r = eye(3) * P1(:,a);
                  if( sqrt((P_r(1) - C1(1))^2 + (P_r(2) - C1(2))^2) < (R1 - w/2) ) && ( (sqrt((P_r(1) - C2(1))^2 + (P_r(2) - C2(2))^2) > (R2 + w/2) || sqrt((P_r(1) - C2(1))^2 + (P_r(2) - C2(2))^2) < (R2 - w/2) || P_r(3) < (C2(3)-w/2) || P_r(3) > (C2(3)+w/2) ))
                    B21(:,a) = integral(@(theta)intgral_func(theta,R2,C2,P1(:,a),I,r),0,2*pi,'ArrayValued',true);
                  end
               end
           end
          M21(ii,di) = sum(B21(3,:)) * dx1 * dx1 / I;
    end
end
% K = M21 / sqrt(L1*L2);
time = toc;
beep
K = M21 / sqrt(L1*L2);
%% plot 
for i=1:N_tilt
    figure('Name','mutual coupling'+string(i))
    plot(RangeD,squeeze(K(i,:))'*100);
    grid on 
    hold on 
    xlabel('X')
    ylabel('Mutual coupling %')
end