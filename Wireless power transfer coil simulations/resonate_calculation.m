clear
Lr = 1.59e-6;%reciver coil mesured (rectangle side 3cm)
Qr = 10;%quality factor mesured at 200k
fmesu = 200e3;
Rr = 2*pi*fmesu*Lr/Qr;
Fres1 = 168/(2*63)*1e6;
C1 = 1/((2*pi*Fres1)^2*Lr)*1e9;

Fres2 = 168/(2*61)*1e6;
C2 = 1/((2*pi*Fres2)^2*Lr)*1e9;

Fres3 = 168/(2*64)*1e6;
C3 = 1/((2*pi*Fres3)^2*Lr)*1e9;
deltaC = C2 - C3;