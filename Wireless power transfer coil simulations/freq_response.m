clear
kt = -0.05;
Rs = 1;
Ls = 4.5e-6;
Cs = 8.2e-9;
M = kt*Ls;
Vs = 1;
phi = -180:5:180;
freq = (700e3:1e3:1e6);
W = 2 * pi * freq;
Z1 = zeros( length(phi) , length(freq) );
Z2 = zeros( length(phi) , length(freq) );
for phi_i = phi
    phi_rad = phi_i * pi / 180;
    V1 = Vs ;
    V2 = Vs * exp(1i*phi_rad);
    Zs = Rs + 1 ./ (1i*W*Cs ) + 1i * W * Ls;
    I1 = Vs * (Zs - exp(1i*phi_rad)*1i * W * M) ./ ( Zs.^2 + (W * M).^2 );
    I2 = Vs * (exp(1i*phi_rad)*Zs - 1i * W * M) ./ ( Zs.^2 + (W * M).^2 );
    % figure
    % plot(freq,abs(I1));
    % hold on 
    % plot(freq,abs(I2));
    Z1(phi_i == phi , :) = V1 ./ I1;
    Z2(phi_i == phi , :) = V2 ./ I2;
end
%%
figure();
subplot(2,2,1)
surf( freq / 1000 , phi ,abs(Z1));
title('Mag Z1');
xlabel('Frequency(KHz)');
ylabel('phase shif(degrees)')
zlim([0 10])
subplot(2,2,3)
surf( freq / 1000 , phi ,angle(Z1)*180/pi);
title('angle Z1');
xlabel('Frequency(KHz)');
ylabel('phase shif(degrees)')
subplot(2,2,2)
surf( freq / 1000 , phi ,abs(Z2));
title('Mag Z2');
xlabel('Frequency(KHz)');
ylabel('phase shif(degrees)')
zlim([0 10])
subplot(2,2,4)
surf( freq / 1000, phi ,angle(Z2)*180/pi);
title('angle Z2');
xlabel('Frequency(KHz)');
ylabel('phase shif(degrees)')
%%
test = 90;
figure('Name',"phi = "+ string(test));
subplot(2,2,1)
plot(freq,abs(Z1(phi == test ,:)));
title('Mag Z1');
subplot(2,2,3)
plot(freq,angle(Z1(phi == test ,:))*180/pi);
title('angle Z1');
subplot(2,2,2)
plot(freq,abs(Z2(phi == test ,:)));
title('Mag Z2');
subplot(2,2,4)
plot(freq,angle(Z2(phi == test ,:))*180/pi);
title('angle Z2');