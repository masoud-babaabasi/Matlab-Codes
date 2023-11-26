function r = intgral_func(theta,R,Center,P,I,Rotation_mat)
u0 = 1.256637e-6;
% dr = [R*cos(theta)+Center(1)-P(1); R*sin(theta)+Center(2)-P(2);Center(3)-P(3)];
dr = [R*cos(theta); R*sin(theta);0];
dr = Rotation_mat * dr;
dr = dr + Center - P;
dl = I*R*[-sin(theta);cos(theta);0];
dl = Rotation_mat * dl;
r = u0 / 4 / pi * cross(dr,dl)/sqrt( dr(1)^2 + dr(2)^2 + dr(3)^2 )^3;
end

