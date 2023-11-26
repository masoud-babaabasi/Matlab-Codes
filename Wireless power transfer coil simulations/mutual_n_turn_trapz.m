clear 
R1 = gpuArray(50e-3)  ; %Meter
turn1 = gpuArray(12);
w = gpuArray(1e-3);
s = gpuArray(0.35e-3);
I = gpuArray(1);
C1 = [ 0 ; 0 ; 0 ];
n_poly1 = gpuArray(6);
XC1 = gpuArray(return_poly_spiral_corners(n_poly1,turn1,R1,w,s));
dx1 = gpuArray(0.01 * R1) ;% dS
Range1 = gpuArray(-(R1/2-w/2):dx1:(R1/2-w/2));
N1 = gpuArray(length(Range1));
P1 = gpuArray([repmat(Range1,1,N1) ; kron(Range1,ones(1,N1)) ; C1(3)*ones(1,N1*N1)]);
B1 = gpuArray(zeros(size(P1)));

R2 = 30e-3 * sqrt(2) ; %Meter
turn2 = 5;
n_poly2 = 4;
XC2 = return_poly_spiral_corners(n_poly2,turn2,R2,w,s);
C2 = [ 0; 0 ; 0 ];
dx2 = 0.05 * R2;% dS = 0.8%
Range2 = -(R2/2-w/2):dx2:(R2/2-w/2);
N2 = length(Range2);
P2 = ([repmat(Range2,1,N2) ; kron(Range2,ones(1,N2)) ; C2(3)*ones(1,N2*N2)]);
B2 = (zeros(size(P2)));
t = linspace(0,1,10000);
%% calculate B1
tic;
L1 = gpuArray(0);
ii = gpuArray(0);
jj = gpuArray(0);
i = gpuArray(0);
j = gpuArray(0);
for ii=1:turn1
  B1 = gpuArray(zeros(size(P1)));
  for j =1:N1
     for i=1:N1
         for jj=1:n_poly1
            a = gpuArray((j-1) + (i-1)*N1 + 1);
            rage_Xc = gpuArray((ii-1)*n_poly1+ 1 : (ii-1)*n_poly1 + n_poly1  + 1); 
            if( is_on_poly(P1(:,a),XC1(rage_Xc,:),w) == 0  && is_inside_poly(P1(:,a),XC1,n_poly1) == 1 ) %not on the poly
                y = intgral_func2(t,XC1((ii-1) * n_poly1 +  jj ,:),XC1((ii-1) * n_poly1 +  jj + 1,:),P1(:,a),I);
                B1(:,a) = B1(:,a) + trapz(t,y,2);
%                 B1(:,a) = B1(:,a) + integral(@(t)intgral_func2(t,XC1((ii-1) * n_poly1 +  jj ,:),XC1((ii-1) * n_poly1 +  jj + 1,:),P1(:,a),I),0,1,'ArrayValued',true);
            end
          end
      end
  end
  
  for jj=1:turn1
      for j =1:N1
          for i=1:N1
               a = (j-1) + (i-1)*N1 + 1;
               rage_Xc = (jj-1)*n_poly1+ 1 : (jj-1)*n_poly1 + n_poly1  + 1; 
               if( is_on_poly(P1(:,a),XC1(rage_Xc,:),w) == 0  && is_inside_poly(P1(:,a),XC1(rage_Xc,:),n_poly1) == 1 ) %not on the poly
                    L1 = L1 + B1(3,a);
               end
          end          
      end
  end
end
u0 = 1.256637e-6;
L1 = L1 * u0 / 4 / pi;
L1 = L1 * dx1 * dx1 / I;
time = toc;