clear
kt = -0.05;
Rs = 1;
Ls = 4.5e-6;
Cs = 8.2e-9;
M = kt*Ls;
phi = 90;
phi_rad = phi * pi / 180;
Vs = 1;
V1 = Vs ;
V2 = Vs * exp(1i*phi_rad);
fres = 1 / ( 2 * pi *sqrt( Ls * Cs) );
fres_in = 1 / ( 2 * pi *sqrt( ( Ls + M ) * Cs) );
fres_anti = 1 / ( 2 * pi *sqrt( ( Ls - M ) * Cs) );
W_c = 2 * pi * fres_in;
Zs = Rs + 1 / (1i*W_c*Cs ) + 1i * W_c * Ls;
I1 = Vs * (Zs - exp(1i*phi_rad)*1i * W_c * M) / ( Zs^2 + (W_c * M)^2 );
I2 = Vs * (exp(1i*phi_rad)*Zs - 1i * W_c * M) / ( Zs^2 + (W_c * M)^2 );
theta1 = (angle(V1) - angle(I1))*180/pi;
theta2 = (angle(V2) - angle(I2))*180/pi;
a1 = atan(( -(1-W_c^2*Ls*Cs)/(W_c*Cs)  - cos(phi_rad)*W_c*M) / (Rs+sin(phi_rad)*W_c*M));
if a1 < 0 && (Rs+sin(phi_rad)*W_c*M) > 0
    a1 = a1 + 2*pi;
elseif a1 < 0 && (Rs+sin(phi_rad)*W_c*M) < 0
    a1 = a1 + pi;
elseif a1 > 0 && (Rs+sin(phi_rad)*W_c*M) < 0
    a1 = a1 + pi;
end
a2 = atan(( -cos(phi_rad)*(1-W_c^2*Ls*Cs)/(W_c*Cs)-W_c*M +Rs*sin(phi_rad)) /(Rs*cos(phi_rad)+sin(phi_rad)*(1-W_c^2*Ls*Cs)/(W_c*Cs)) );
if a2 < 0 && (Rs*cos(phi_rad)+sin(phi_rad)*(1-W_c^2*Ls*Cs)/(W_c*Cs)) > 0
    a2 = a2 + 2*pi;
elseif a2 < 0 && (Rs*cos(phi_rad)+sin(phi_rad)*(1-W_c^2*Ls*Cs)/(W_c*Cs)) < 0
    a2 = a2 + pi;
elseif a1 > 0 && (Rs*cos(phi_rad)+sin(phi_rad)*(1-W_c^2*Ls*Cs)/(W_c*Cs)) < 0
    a2 = a2 + pi;
end
t1 =(a1 - a2)*180/pi;
t2 = (angle(I1) - angle(I2))*180/pi ;
%%
R = 25e-3 ; %Meter
w = 1e-3;
s = 0.35e-3;

dx = R/3;
err = 10e-3;

w = w / R;
s = s / R;
dx = dx / R;
err = err / R;
R = 1;
Range = -3*R:dx:3*R;
N = length(Range);
P = [zeros(1,N*N) ; repmat(Range,1,N) ; kron(Range,ones(1,N))];

B1 = zeros(size(P));
B2 = zeros(size(P));
Bt = zeros([length(phi_rad) , size(P)]);
C1 = [ 0 ; -R ; 0 ];
C2 = [ 0 ; R ; 0 ];
for i_phase = 1:length(phi_rad)
    for i=1:N
       for j =1:N
              a =  (j-1) + (i-1)*N + 1;
              if (P(1,a)^2 + (P(2,a)-C1(2))^2 + P(3,a)^2 + R^2 - 2*R*sqrt(P(1,a)^2 + (P(2,a)-C1(2))^2)) > (err)^2 && (P(1,a)^2 + (P(2,a)-C2(2))^2 + P(3,a)^2 + R^2 - 2*R*sqrt(P(1,a)^2 + (P(2,a)-C2(2))^2)) > (err)^2 
                  B1(:,a) = integral(@(theta)intgral_func(theta,R,C1,P(:,a),I1,eye(3)),0,2*pi,'ArrayValued',true);
                  B2(:,a) = integral(@(theta)intgral_func(theta,R,C2,P(:,a),I2,eye(3)),0,2*pi,'ArrayValued',true);
                  Bt(i_phase,:,:) = B1 + B2;
              end
       end
    end
end
%%
for i_phase = 1:length(phi_rad)
    figure('Name','field')
    quiver(P(2,:),P(3,:),abs(squeeze(Bt(i_phase,2,:)))'.*sign(cos(angle(squeeze(Bt(i_phase,2,:)))))',abs(squeeze(Bt(i_phase,3,:)))'.*sign(cos(angle(squeeze(Bt(i_phase,3,:)))))','B')
    hold on 
    marker_size = 10;
    plot(2*R,0,'o','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)
    plot(2*R,0,'.','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)
    plot(w+s,0,'o','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)
    plot(w+s,0,'+','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)

    plot(-(w+s),0,'o','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
    plot(-(w+s),0,'+','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
    plot(-2*R,0,'o','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
    plot(-2*R,0,'.','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
    xlim([-3 3])
    ylim([-2 2])
    xlabel('Y')
    ylabel('Z')
    title("Phase shift = " + string(phi_rad(i_phase)*180/pi));
end
%%
figure('Name','field')
hold on
plot(2*R,0,'o','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(2*R,0,'.','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(w+s,0,'o','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)
plot(w+s,0,'+','MarkerFaceColor','white','MarkerEdgeColor','blue','MarkerSize',marker_size)

plot(-(w+s),0,'o','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
plot(-(w+s),0,'+','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
plot(-2*R,0,'o','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
plot(-2*R,0,'.','MarkerFaceColor','white','MarkerEdgeColor','red','MarkerSize',marker_size)
xlim([-3 3])
ylim([-2 2])
xlabel('Y')
ylabel('Z')
marker_size = 10;
for i_phase = 1:length(phi_rad)
    quiver(P(2,219),P(3,219),abs(squeeze(Bt(i_phase,2,219)))*sign(cos(angle(squeeze(Bt(i_phase,2,219)))))*0.25e7,abs(squeeze(Bt(i_phase,3,219)))*sign(cos(angle(squeeze(Bt(i_phase,3,219)))))*0.25e7,'B')
    hold on 
end