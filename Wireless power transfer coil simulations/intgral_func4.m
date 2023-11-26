function r = intgral_func4(theta,R,Center,w,s,l,P,I)
u0 = 1.256637e-6;
dr = [R*cos(theta)+Center(1)-P(1); R*sin(theta)+Center(2)-P(2);-l/2+theta/2/pi*(w+s)+Center(3)-P(3)];
dl = I*R*[-sin(theta);cos(theta);(w+s)/2/pi/R];
r = u0 / 4 / pi * cross(dr,dl)/sqrt( dr(1)^2 + dr(2)^2 + dr(3)^2 )^3;
end

