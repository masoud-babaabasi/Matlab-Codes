clear
dout = 30e-3;
% din = (30:80)* 1e-3;
din = 11e-3;
f = 10e3:1e3:1e6;
% f = 1e6;
s = 350e-6;
w = 1e-3;

c1 = 1.27;
c2 = 2.07;
c3 = 0.18;
c4 = 0.13;

u = 4*pi*1e-7 ;
e = 2.718281;

t= 18e-6;
r0 = 1.72e-8;

davg = (dout+din)/2;
ro = (dout - din)/(dout + din);

n = (dout - din + 2*s)/(2*(w+s));

% lres = dout + 4*n*dout - 2*n.*(2*n-1)*(w+s)
lres = Spiral_length(dout,w,s,n,length(n));


L = 1/2*u*n.^2.*davg*c1*(log(c2/ro)/log(e)+c3*ro+c4*ro^2);
res = r0 * lres / (t * w);
% not included the stray capacitance
Q = 2*pi*f * L ./ res ;

% figure('Name','Indutance','NumberTitle','off'); 
% plot(din,L)
% figure('Name','Quality factor','NumberTitle','off'); 
% plot(din,Q)
% figure('Name','Number of turns','NumberTitle','off'); 
% plot(din,n)

figure('Name','Quality factor','NumberTitle','off'); 
plot(f,Q)