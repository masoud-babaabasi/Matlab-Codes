rmean = 0;
r1 = rmean + 6* sin(2*pi/16*(1:32));
r2 = rmean + 6* sin(2*pi/16*(1:2^10));
rf1 = fft(r1 / length(r1));
f1 = 0:1/length(r1):1-1/length(r1);
rf2 = fft(r2 / length(r2));
f2 = 0:1/length(r2):1-1/length(r2);
plot( f1(1:length(r1)/2) , abs(rf1(1:length(r1)/2)) )
hold on 
plot( f2(1:length(r2)/2) , abs(rf2(1:length(r2)/2)) ,'r')
xlabel('nirmalazed frequency')