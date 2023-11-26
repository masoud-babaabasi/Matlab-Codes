clear 
R = 25e-3 ; %Meter
I = 1;
P = [0 ; 0 ; 30]*1e-3;
C1 = [ 0 ; -R ; 0 ];
C2 = [ 0 ; R; 0 ];
%**********************
phi = 0:pi/8:2*pi;
n = length(phi);
B = zeros(3,n);
for i = 1:n
    if sin(phi(i)) ~= 0
        q1 = integral(@(theta)intgral_func(theta,R,C1,P,I*sin(phi(i))),0,2*pi,'ArrayValued',true);
    else 
        q1 = [0;0;0];
    end
    if cos(phi(i)) ~= 0
        q2 = integral(@(theta)intgral_func(theta,R,C2,P,I*cos(phi(i))),0,2*pi,'ArrayValued',true);
    else
        q2 = [0;0;0];
    end
    B(:,i) = q1 + q2; 
end
B = B ./ vecnorm(B) * R;
%%
t = 0:1e-3:2*pi;
coil1 = [R*cos(t)+C1(1); R*sin(t)+C1(2);ones(1,length(t))*(C1(3))];
coil2 = [R*cos(t)+C2(1); R*sin(t)+C2(2);ones(1,length(t))*(C2(3))];
plot3(coil1(1,:),coil1(2,:),coil1(3,:))
hold on
plot3(coil2(1,:),coil2(2,:),coil2(3,:))
quiver3(ones(1,n)*P(1),ones(1,n)*P(2),ones(1,n)*P(3),B(1,:),B(2,:),B(3,:))
axis equal