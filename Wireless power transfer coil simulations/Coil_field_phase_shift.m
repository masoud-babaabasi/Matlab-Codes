clear 
R = 25e-3 ; %Meter
w = 1e-3;
s = 0.35e-3;
I = 1;
% dx = 10e-3;
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
phase_shift = 0:pi/6:pi;
%%
B1 = zeros(size(P));
B2 = zeros(size(P));
Bt = zeros([length(phase_shift) , size(P)]);
C1 = [ 0 ; -R ; 0 ];
C2 = [ 0 ; R ; 0 ];
for i_phase = 1:length(phase_shift)
    for i=1:N
       for j =1:N
              a =  (j-1) + (i-1)*N + 1;
              if (P(1,a)^2 + (P(2,a)-C1(2))^2 + P(3,a)^2 + R^2 - 2*R*sqrt(P(1,a)^2 + (P(2,a)-C1(2))^2)) > (err)^2 && (P(1,a)^2 + (P(2,a)-C2(2))^2 + P(3,a)^2 + R^2 - 2*R*sqrt(P(1,a)^2 + (P(2,a)-C2(2))^2)) > (err)^2 
                  B1(:,a) = integral(@(theta)intgral_func(theta,R,C1,P(:,a),I,eye(3)),0,2*pi,'ArrayValued',true);
                  B2(:,a) = integral(@(theta)intgral_func(theta,R,C2,P(:,a),I*exp(phase_shift(i_phase)*sqrt(-1)),eye(3)),0,2*pi,'ArrayValued',true);
                  Bt(i_phase,:,:) = B1 + B2;
              end
       end
    end
end
%%
for i_phase = 1:length(phase_shift)
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
    title("Phase shift = " + string(phase_shift(i_phase)*180/pi));
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
for i_phase = 1:length(phase_shift)
    quiver(P(2,219),P(3,219),abs(squeeze(Bt(i_phase,2,219)))*sign(cos(angle(squeeze(Bt(i_phase,2,219)))))*0.25e7,abs(squeeze(Bt(i_phase,3,219)))*sign(cos(angle(squeeze(Bt(i_phase,3,219)))))*0.25e7,'B')
    hold on 
end