function df = intgral_func3(theta,Rin,w,s,C,p,I)
u0 = 1.256637e-6;
R = Rin + theta * ( w + s ) / 2 /pi ;
dl = I*R*[-sin(theta);cos(theta);0];
x = [R*cos(theta)+C(1);R*sin(theta)+C(2);C(3)];
r = [p(1)-x(1) ; p(2)-x(2) ; p(3)-x(3)];
df = u0 / 4 / pi * cross(dl,r)/sqrt( r(1)^2 + r(2)^2 + r(3)^2 )^3;
end