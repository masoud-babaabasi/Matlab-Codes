clear 
R = 25e-3 ; %Meter
w = 1e-3;
I = 1;
C1 = [ 0 ; 0 ; 0 ];
C2 = [ 1e-3; 0 ; 0 ];
dx = 0.01 * R;% dS = 0.8%
Range = -(R-w/2):dx:(R-w/2);
N = length(Range);
P = [repmat(Range,1,N) ; kron(Range,ones(1,N)) ; zeros(1,N*N)];
B1 = zeros(size(P));
B2 = zeros(size(P));

%%
tic;
for i=1:N
   for j =1:N
      a = (j-1) + (i-1)*N + 1;
      if( sqrt((P(1,a) - C1(1))^2 + (P(2,a) - C1(2))^2) < (R - w/2) )
        B1(:,a) = integral(@(theta)intgral_func(theta,R,C1,P(:,a),I,eye(3)),0,2*pi,'ArrayValued',true);
      end
   end
end
L = abs(sum(B1(3,:))) * dx * dx / I;
%%
D =0:R/15:3*R;
M21 = zeros(1,length(D));
for d = D
   C2 = [ 0; d ; 0];
   for i=1:N
       for j =1:N
          a = (j-1) + (i-1)*N + 1;
          if( sqrt((P(1,a) - C1(1))^2 + (P(2,a) - C1(2))^2) < (R - w/2)  && ( sqrt((P(1,a) - C2(1))^2 + (P(2,a) - C2(2))^2) > (R + w/2) || sqrt((P(1,a) - C2(1))^2 + (P(2,a) - C2(2))^2) < (R - w/2)) )
            B2(:,a) = integral(@(theta)intgral_func(theta,R,C2,P(:,a),I,eye(3)),0,2*pi,'ArrayValued',true);
            
%             if( abs(B2(3,a)) >= max(abs(B1(3,:))) )
%                 %phi12(1,find(D==d)) = phi12(1,find(D==d)) + B2(3,a) * dx * dx ;
%                 B2(3,a) = 0;
%             end
          end
       end
   end
  M21(1,find(D==d)) = sum(B2(3,:)) * dx * dx / I;
%     disp(phi12(1,find(D==d)))
disp(d)
end

time = toc;
%%
K = M21 / L;
figure('Name','One coil')
quiver3(P(1,:),P(2,:),P(3,:),B1(1,:),B1(2,:),B1(3,:),'R')
hold on
B2(1,abs(B2(1,:)) > max(abs(B1(3,:))) ) = 0;
B2(2,abs(B2(2,:)) > max(abs(B1(3,:))) ) = 0;
B2(3,abs(B2(3,:)) > max(abs(B1(3,:))) ) = 0;
quiver3(P(1,:),P(2,:),P(3,:),B2(1,:),B2(2,:),B2(3,:),'B')

t = 0:1e-3:2*pi;
coil1 = [R*cos(t)+C1(1); R*sin(t)+C1(2);ones(1,length(t))*(C1(3))];
coil2 = [R*cos(t)+C2(1); R*sin(t)+C2(2);ones(1,length(t))*(C2(3))];
plot3(coil1(1,:),coil1(2,:),coil1(3,:))
plot3(coil2(1,:),coil2(2,:),coil2(3,:))
dl = I*[-sin(t(1));cos(t(1));0]*R;
quiver3(coil1(1),coil1(2),coil1(3),dl(1),dl(2),dl(3),'G')
figure('Name','mutual coupling')
plot(D/R,K)
grid on 
xlim([0 3])
xlabel('Distance')
xticklabels([ "0" string(0.5:0.5:3)+'R'])
ylabel('Mutual coupling')
ax = gca;
grid on
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';