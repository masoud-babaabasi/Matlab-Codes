clear 
R1 = 50e-3/2 ; %Meter
w = 1e-3;
I = 1;
C1 = [ 0 ; 0 ; 0 ];
dx1 = 0.02 * R1;% dS = 0.8%
Range1 = -(R1-w/2):dx1:(R1-w/2);
N1 = length(Range1);
P1 = [repmat(Range1,1,N1) ; kron(Range1,ones(1,N1)) ; C1(3)*ones(1,N1*N1)];
B1 = zeros(size(P1));

R2 = 30e-3/2 ; %Meter
C2 = [ 0; 0 ; 40e-3 ];
dx2 = 0.02 * R2;% dS = 0.8%
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
Range_tilt = 0:pi/18:pi/2;
N_tilt = length(Range_tilt);
M21 = zeros(N_tilt,3);
P = [0.5*R1 , -R1 ; 0.5*R1 , R1 ;-(sqrt(3) - 0.5)*R1 , 0 ];
for ii = 1:N_tilt
    x_angle = 0;
    y_angle = Range_tilt(ii);
    z_angle = 0;
    r_z = [cos(z_angle) -sin(z_angle) 0 ; sin(z_angle) cos(z_angle) 0 ;0 0 1];
    r_y = [cos(y_angle) 0 sin(y_angle)  ; 0 1 0 ; -sin(y_angle) 0 cos(y_angle)];
    r_x = [1 0 0; 0 cos(x_angle) -sin(x_angle);0 sin(x_angle) cos(x_angle) ];
    r = r_x*r_y*r_z;
    for di = 1:3
       C2 = [ P(di,1); P(di,2) ; 40e-3];
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
       %B21 = r*B21;
      M21(ii,di) = sum(B21(3,:)) * dx1 * dx1 / I;
    end
end
% K = M21 / sqrt(L1*L2);
time = toc;
beep
K = M21 / sqrt(L1*L2);
%% plot
for i=1:3
    figure
    plot(K(:,i))
    grid on 
    hold on 
    xlabel('X')
    ylabel('Mutual coupling %')
    % xlim([-3 3])
    % xticks([-3:1:3])
end
%%
n = sqrt(L2/L1);
Rs = 1;
RL = 500;
Vs = 1;
power1 = (K(:,1)/n).^2*RL ./ (Rs + (K(:,1)/n).^2*RL).^2 * Vs^2;
power2 = (K(:,2)/n).^2*RL ./ (Rs + (K(:,2)/n).^2*RL).^2 * Vs^2;
power3 = (K(:,3)/n).^2*RL ./ (Rs + (K(:,3)/n).^2*RL).^2 * Vs^2;
power_in = ( sqrt(power1) + sqrt(power2) + sign(K(:,3)).*sqrt(power3)).^2;
power_an = ( sqrt(power1) + sqrt(power2) - sign(K(:,3)).*sqrt(power3)).^2;

effi1 = (K(:,1)/n).^2*RL ./ (Rs + (K(:,1)/n).^2*RL);
effi2 = (K(:,2)/n).^2*RL ./ (Rs + (K(:,2)/n).^2*RL);
effi3 = (K(:,3)/n).^2*RL ./ (Rs + (K(:,3)/n).^2*RL);
effi_anti_mid = ( sqrt(effi1) + sqrt(effi2) - sign(K(:,3)).*sqrt(effi3)).^2/3;
effi_in_mid = ( sqrt(effi1) + sqrt(effi2) + sign(K(:,3)).*sqrt(effi3)).^2/3;
figure
plot(Range_tilt * 180 / pi , effi_in_mid.*power_in )
hold on 
plot(Range_tilt * 180 / pi , effi_anti_mid.*power_an )
xlabel('Receiver Coil Angle( \circ )')
ylabel('FOM')
legend(["Inphase" "Antiphase"])
grid on