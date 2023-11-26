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
D = -1*R1:2*R1:1*R1;
Range_tilt = 0:pi/36:pi/2;
N_tilt = length(Range_tilt);
tilt = [repmat(Range_tilt,1,N_tilt) ; kron(Range_tilt,ones(1,N_tilt))];
M21 = zeros(N_tilt,length(D));
for ii = 1:N_tilt
        x_angle = Range_tilt(ii);
        y_angle = 0;
        z_angle = 0;
        r_z = [cos(z_angle) -sin(z_angle) 0 ; sin(z_angle) cos(z_angle) 0 ;0 0 1];
        r_y = [cos(y_angle) 0 sin(y_angle)  ; 0 1 0 ; -sin(y_angle) 0 cos(y_angle)];
        r_x = [1 0 0; 0 cos(x_angle) -sin(x_angle);0 sin(x_angle) cos(x_angle) ];
        r = r_x*r_y*r_z;
        for d = D
           C2 = [ 0; d ; 40e-3];
           B21 = zeros(size(P1));
           for i=1:N1
               for j =1:N1
                  a = (j-1) + (i-1)*N1 + 1;
                  P_r = eye(3) * P1(:,a);
                  if( sqrt((P_r(1) - C1(1))^2 + (P_r(2) - C1(2))^2) < (R1 - w/2) ) && ( (sqrt((P_r(1) - C2(1))^2 + (P_r(2) - C2(2))^2) > (R2 + w/2) || sqrt((P_r(1) - C2(1))^2 + (P_r(2) - C2(2))^2) < (R2 - w/2) || P_r(3) < (C2(3)-w/2) || P_r(3) > (C2(3)+w/2) ))
                    B21(:,a) = integral(@(theta)intgral_func(theta,R2,C2,P1(:,a),I,r),0,2*pi,'ArrayValued',true);

        %             if( abs(B2(3,a)) <= max(abs(B1(3,:))) )
        %                 phi12(1,find(D==d)) = phi12(1,find(D==d)) + B2(3,a) * dx * dx ;
        %             end
                  end
               end
           end
           %B21 = r*B21;
          M21(ii,D==d) = sum(B21(3,:)) * dx1 * dx1 / I;
        %     disp(phi12(1,find(D==d)))
        end
end
% K = M21 / sqrt(L1*L2);
time = toc;
beep
K = M21 / sqrt(L1*L2);
%%
% figure('Name','mutual coupling')
% plot(D/R1,K'*100)
% grid on 
% legend('x='+string(round(tilt(1,:)*180/pi))+'y='+string(round(tilt(2,:)*180/pi)))
% xlabel('Normalized distance (D/R)')
% ylabel('Mutual coupling %')
%% plot 
for i=2:2
    figure('Name','mutual coupling'+string(i))
    plot(D/R1,K(i,:)*100)
    grid on 
    hold on 
    plot(D/R1+2,K(i,:)*100)
%     plot(D(D/R1 >= -3 & D/R1<=5)/R1,( K(1,D/R1 >= -3 & D/R1<=5) + K(1,D/R1 >= -5 & D/R1<=3))*100)
%     plot(D/R1-2,K(i,:)*100)
    legend(['K1r';'K2r'])
    xlim([-3 5])
    %legend('x='+string(round(tilt(1,i:N_tilt:i+N_tilt*(N_tilt-1))*180/pi))+' ,y='+string(round(tilt(2,i:N_tilt:i+N_tilt*(N_tilt-1))*180/pi)))
    xlabel('Normalized distance (D/R)')
    ylabel('Mutual coupling %')
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
end
%% power and effi for parallel coil
n = sqrt(L2/L1);
Rs = 1;
RL = 100;
Vs = 1;
effi = (K(1,:)/n).^2*RL ./ (Rs + (K(1,:)/n).^2*RL);
power = (K(1,:)/n).^2*RL ./ (Rs + (K(1,:)/n).^2*RL).^2 * Vs^2;
figure
plot(D/R1,effi*100)
hold on 
plot(D/R1+2,effi*100)
plot(D(D/R1 >= -2 & D/R1<=4)/R1,( effi(D/R1 >= -2 & D/R1<=4) + effi(D/R1 >= -4 & D/R1<=2))*100/2)
figure
plot(D/R1,power/max(power)*100)
hold on
plot(D/R1+2,power/max(power)*100)
power_inphase = ( sqrt(power(D/R1 >= -3 & D/R1<=5)) + sqrt(power(D/R1 >= -5 & D/R1<=3))).^2;
plot(D(D/R1 >= -3 & D/R1<=5)/R1,power_inphase*100/max(power))
%% FOM
figure
plot(D/R1,power.*effi)
hold on 
plot(D/R1+2,power.*effi)
effi_inphase = ( sqrt(effi(D/R1 >= -3 & D/R1<=5)) + sqrt(effi(D/R1 >= -5 & D/R1<=3))).^2/2;
plot(D(D/R1 >= -3 & D/R1<=5)/R1,power_inphase.*effi_inphase)
legend(["Coil1"  "Coil2" "Both coils"])
xlim([-2 4])
xlabel('Normalized distance (D/R)')
ylabel('FOM')
ax = gca;
xticks(-2:0.5:4)
grid on
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
%%
% effi2 = (K(:,41)/n).^2*RL ./ (Rs + (K(:,41)/n).^2*RL);
power1 = (K(:,D==R1)/n).^2*RL ./ (Rs + (K(:,D==R1)/n).^2*RL).^2 * Vs^2;
power2 = (K(:,D==-R1)/n).^2*RL ./ (Rs + (K(:,D==-R1)/n).^2*RL).^2 * Vs^2;
effi1 = (K(:,D==R1)/n).^2*RL ./ (Rs + (K(:,D==R1)/n).^2*RL);
effi2 = (K(:,D==-R1)/n).^2*RL ./ (Rs + (K(:,D==-R1)/n).^2*RL);
% figure
% plot(Range_tilt*180/pi,effi2*100)
% figure
% plot(Range_tilt*180/pi,power2/max(power2)*100)
power_in_mid = ( sqrt(power2) + sign(K(:,D==R1)).*sqrt(power1)).^2;
power_anti_mid = ( sqrt(power2) - sign(K(:,D==R1)).*sqrt(power1)).^2;
figure
plot(Range_tilt*180/pi,power_in_mid)
hold on 
plot(Range_tilt*180/pi,power_anti_mid)
effi_anti_mid = ( sqrt(effi2) - sign(K(:,D==R1)).*sqrt(effi1)).^2/2;
effi_in_mid = ( sqrt(effi2) + sign(K(:,D==R1)).*sqrt(effi1)).^2/2;
figure
plot(Range_tilt*180/pi,effi_in_mid)
hold on 
plot(Range_tilt*180/pi,effi_anti_mid)
%% FOM
figure
plot(Range_tilt*180/pi,effi_in_mid.*power_in_mid)
hold on 
plot(Range_tilt*180/pi,effi_anti_mid.*power_anti_mid)
xlabel('Receiver Coil Angle( \circ )')
ylabel('FOM')
legend(["Inphase" "Antiphase"])
ax = gca;
grid on
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
%%
% surf(D/R1,Range_tilt*180/pi,K*100,'FaceColor','interp');
% xlabel('Normalized distance (D/R)')
% ylabel('Angle')
% ylim([0 90])
% zlabel('Mutual coupling %')