function df = intgral_func2(t,x1,x2,p,I)
u0 = 1.256637e-6;
dl = [x2(1)-x1(1);x2(2)-x1(2);x2(3)-x1(3)];
x = [dl(1)*t+x1(1);dl(2)*t+x1(2);dl(3)*t+x1(3)];
r = [p(1)-x(1,:) ; p(2)-x(2,:) ; p(3)-x(3,:)];
%dl = repmat(dl,1,length(x));
df = u0 / 4 / pi * I * cross(dl,r)/sqrt( r(1)^2 + r(2)^2 + r(3)^2 )^3;
end
