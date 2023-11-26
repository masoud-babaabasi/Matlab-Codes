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
for i=1:N2
   for j =1:N2
      a = (j-1) + (i-1)*N2 + 1;
      if( sqrt((P2(1,a) - C2(1))^2 + (P2(2,a) - C2(2))^2) < (R2 - w/2) )
        B2(:,a) = integral(@(theta)intgral_func(theta,R2,C2,P2(:,a),I,eye(3)),0,2*pi,'ArrayValued',true);
      end
   end
end
L2 = abs(sum(B2(3,:))) * dx2 * dx2 / I;
rangeR = R1:R1/5:8*R1;
r = eye(3);
M21 = zeros(1,length(rangeR));
K = zeros(1,length(rangeR));
for R1 = rangeR  %Meter
    dx1 = 0.03 * R1;% dS = 0.8%
    Range1 = -(R1-w/2):dx1:(R1-w/2);
    N1 = length(Range1);
    P1 = [repmat(Range1,1,N1) ; kron(Range1,ones(1,N1)) ; C1(3)*ones(1,N1*N1)];
    B1 = zeros(size(P1));
    for i=1:N1
       for j =1:N1
          a = (j-1) + (i-1)*N1 + 1;
          if( sqrt((P1(1,a) - C1(1))^2 + (P1(2,a) - C1(2))^2) < (R1 - w/2) )
            B1(:,a) = integral(@(theta)intgral_func(theta,R1,C1,P1(:,a),I,eye(3)),0,2*pi,'ArrayValued',true);
          end
       end
    end
    L1 = abs(sum(B1(3,:))) * dx1 * dx1 / I;
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
%   M21(R1 == rangeR) = sum(B21(3,:)) * dx1 * dx1 / I;
  K(R1 == rangeR) = sum(B21(3,:)) * dx1 * dx1 / I / sqrt(L1*L2);
end
% K = M21 / sqrt(L1*L2);
time = toc;
beep


%% plot 
    figure('Name','mutual coupling'+string(i))
    plot(rangeR/R2,K*100)
    grid on 
    hold on 
%     plot(D/R1-2,K(i,:)*100)
    %legend('x='+string(round(tilt(1,i:N_tilt:i+N_tilt*(N_tilt-1))*180/pi))+' ,y='+string(round(tilt(2,i:N_tilt:i+N_tilt*(N_tilt-1))*180/pi)))
    xlabel('Tx coil radius (R1/R2)')
    ylabel('Mutual coupling(K12) %')
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';