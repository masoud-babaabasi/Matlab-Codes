clear 
load('mutual3.mat')
clear P1 P2
R1 = 50e-3/2 ; %Meter
R = 1;
C1 = [ - R /2 ; +R ; 0 ];
C2 = [ - R /2 ; -R ; 0 ];
C3 = [1.7*R - R /2 ; 0 ; 0 ];
K = squeeze(K(1,:,:));
%%
D = -3:0.2:3;
figure
surf(D( D >= -1.5 & D <= 2),D(D >= -2 & D <= 2),K((D+C1(2)) >= -2 & (D+C1(2)) <= 2 , (D+C1(1)) >= -1.5 & (D+C1(1)) <= 2)*100,'FaceColor','interp');
xlabel('X')
ylabel('Y')
hold on 
surf(D( D >= -1.5 & D <= 2),D(D >= -2 & D <= 2),K((D+C2(2)) >= -2 & (D+C2(2)) <= 2 , (D+C2(1)) >= -1.5 & (D+C2(1)) <= 2)*100,'FaceColor','interp');
surf(D( D >= -1.5 & D <= 2),D(D >= -2 & D <= 2),K((D+C3(2)) >= -2 & (D+C3(2)) <= 2 , (D+C3(1)) >= -1.5 & (D+C3(1)) <= 2)*100,'FaceColor','interp');
%%
rangeX = D( D >= -1.5 & D <= 2);
rangeY = D(D >= -2 & D <= 2);
K1 = K((D+C1(2)) >= -2 & (D+C1(2)) <= 2 , (D+C1(1)) >= -1.5 & (D+C1(1)) <= 2) ; 
K2 = K((D+C2(2)) >= -2 & (D+C2(2)) <= 2 , (D+C2(1)) >= -1.5 & (D+C2(1)) <= 2) ; 
K3 = K((D+C3(2)) >= -2 & (D+C3(2)) <= 2 , (D+C3(1)) >= -1.5 & (D+C3(1)) <= 2) ; 
figure
surf(rangeX , rangeY , K1);
hold on 
surf(rangeX , rangeY , K2);
surf(rangeX , rangeY , K3);
%% FOM
n = sqrt(L2/L1);
Rs = 1;
RL = 500;
Vs = 1;
power = zeros([ 7,size(K1)]);
P1 = (K1/n).^2*RL ./ (Rs + (K1/n).^2*RL).^2 * Vs^2;
P2 = (K2/n).^2*RL ./ (Rs + (K2/n).^2*RL).^2 * Vs^2;
P3 = (K3/n).^2*RL ./ (Rs + (K3/n).^2*RL).^2 * Vs^2;
power(1,:,:) = P1;
power(2,:,:) = P2;
power(3,:,:) = P3;
power(4,:,:) = ( sqrt(P1).*sign(K1) +  sqrt(P2).*sign(K2) ).^2;
power(5,:,:) = ( sqrt(P1).*sign(K1) +  sqrt(P3).*sign(K3) ).^2;
power(6,:,:) = ( sqrt(P2).*sign(K2) +  sqrt(P3).*sign(K3) ).^2;
power(7,:,:) = ( sqrt(P1).*sign(K1) +  sqrt(P2).*sign(K2) + sqrt(P3).*sign(K3) ).^2;

effi = zeros([ 7,size(K1)]);
e1 = (K1/n).^2*RL ./ (Rs + (K1/n).^2*RL);
e2 = (K2/n).^2*RL ./ (Rs + (K2/n).^2*RL);
e3 = (K3/n).^2*RL ./ (Rs + (K3/n).^2*RL);
effi(1,:,:) = e1;
effi(2,:,:) = e2;
effi(3,:,:) = e3;
effi(4,:,:) = ( sqrt(e1).*sign(K1) +  sqrt(e2).*sign(K2) ).^2/2;
effi(5,:,:) = ( sqrt(e1).*sign(K1) +  sqrt(e3).*sign(K3) ).^2/2;
effi(6,:,:) = ( sqrt(e2).*sign(K2) +  sqrt(e3).*sign(K3) ).^2/2;
effi(7,:,:) = ( sqrt(e1).*sign(K1) +  sqrt(e2).*sign(K2) +  sqrt(e3).*sign(K3) ).^2/3;
FOM = power .* effi;
% figure
% surf(rangeX , rangeY , squeeze(effi(4,:,:)));
% hold on 
% for i=5:6
%     surf(rangeX , rangeY , squeeze(effi(i,:,:)));
% end
marker_size = 8;
figure 
hold on
plot(-5,-5,'x','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(-5,-5,'O','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
plot(-5,-5,'*','MarkerFaceColor','white','MarkerEdgeColor','green','MarkerSize',marker_size)
for i = 1:length(rangeX)
   for j = 1:length(rangeY)
       [f , inx ] =max( effi(:,j , i));
       if  inx <= 3 
           mark = 'x';
           color = 'blue';
       elseif inx <= 6 && inx >= 4
           mark = 'O';
           color = 'red';
       else 
           mark = '*';
           color = 'green';
       end
       plot(rangeX(i),rangeY(j),mark,'MarkerFaceColor','white','MarkerEdgeColor',color,'MarkerSize',marker_size)
   end
end
axis equal
xlabel('X')
ylabel('Y')

xlim([-1.5 2])
ylim([-2 2])

t = 0:1e-3:2*pi;
C3 = [sqrt(3)*R - R /2 ; 0 ; 0 ];
coil1 = [R*cos(t)+C1(1); R*sin(t)+C1(2)];
coil2 = [R*cos(t)+C2(1); R*sin(t)+C2(2)];
coil3 = [R*cos(t)+C3(1); R*sin(t)+C3(2)];
plot(coil1(1,:),coil1(2,:),'black')
plot(coil2(1,:),coil2(2,:),'black')
plot(coil3(1,:),coil3(2,:),'black')
plot(C1(1),C1(2),'.','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',marker_size)
plot(C2(1),C2(2),'.','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',marker_size)
plot(C3(1),C3(2),'.','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',marker_size)
legend(["Single Transmitter" , "Dual transmitter inphase" , "Tripple transmitter inphase"],'Location','northeastoutside')
