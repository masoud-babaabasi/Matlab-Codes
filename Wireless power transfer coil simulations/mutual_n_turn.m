clear 
R1 = 50e-3  ; %Meter
turn1 = 12;
w = 1e-3;
s = 0.35e-3;
I = 1;
C1 = [ 0 ; 0 ; 0 ];
n_poly1 = 6;
XC1 = return_poly_spiral_corners(n_poly1,turn1,R1,w,s);
dx1 = 0.05 * R1 ;% dS
Range1 = (-(R1/2-w/2):dx1:(R1/2-w/2));
N1 = length(Range1);
P1 = ([repmat(Range1,1,N1) ; kron(Range1,ones(1,N1)) ; C1(3)*ones(1,N1*N1)]);
B1 = (zeros(size(P1)));

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
%% calculate L1
tic;
L1 = 0;
for ii=1:turn1
  B1 = zeros(size(P1));
  for j =1:N1
     for i=1:N1
         for jj=1:n_poly1
            a = (j-1) + (i-1)*N1 + 1;
            rage_Xc = (ii-1)*n_poly1+ 1 : (ii-1)*n_poly1 + n_poly1  + 1; 
            if( is_on_poly(P1(:,a),XC1(rage_Xc,:),w) == 0  && is_inside_poly(P1(:,a),XC1,n_poly1) == 1 ) %not on the poly
                B1(:,a) = B1(:,a) + integral(@(t)intgral_func2(t,XC1((ii-1) * n_poly1 +  jj ,:),XC1((ii-1) * n_poly1 +  jj + 1,:),P1(:,a),I),0,1,'ArrayValued',true);
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
L1 = L1 * dx1 * dx1 / I;
%% calculate L2
L2 = 0;
for ii=1:turn2
  B2 = zeros(size(P2));
  for j =1:N2
     for i=1:N2
         for jj=1:n_poly2
            a = (j-1) + (i-1)*N2 + 1;
            rage_Xc = (ii-1)*n_poly2+ 1 : (ii-1)*n_poly2 + n_poly2  + 1; 
            if( is_on_poly(P2(:,a),XC2(rage_Xc,:),w) == 0  && is_inside_poly(P2(:,a),XC2,n_poly2) == 1 ) %not on the poly
                B2(:,a) = B2(:,a) + integral(@(t)intgral_func2(t,XC2((ii-1) * n_poly2 +  jj ,:),XC2((ii-1) * n_poly2 +  jj + 1,:),P2(:,a),I),0,1,'ArrayValued',true);
            end
          end
      end
  end
  
  for jj=1:turn2
      for j =1:N2
          for i=1:N2
               a = (j-1) + (i-1)*N2 + 1;
               rage_Xc = (jj-1)*n_poly2+ 1 : (jj-1)*n_poly2 + n_poly2  + 1; 
               if( is_on_poly(P2(:,a),XC2(rage_Xc,:),w) == 0  && is_inside_poly(P2(:,a),XC2(rage_Xc,:),n_poly2) == 1 ) %not on the poly
                    L2 = L2 + B2(3,a);
               end
          end          
      end
  end
end
L2 = L2 * dx2 * dx2 / I;
%% claculate K
RangeX = -4*R1:5e-3:4*R1;
Nx = length(RangeX);
Range_tilt = 0:pi/18:pi/2;
N_tilt = length(Range_tilt);

M12 = zeros(N_tilt,Nx,Nx);
parfor i_t = 1:N_tilt
    Range = RangeX;
    x_angle = Range_tilt(i_t);
    y_angle = 0;
    z_angle = 0;
    r_z = [cos(z_angle) -sin(z_angle) 0 ; sin(z_angle) cos(z_angle) 0 ;0 0 1];
    r_y = [cos(y_angle) 0 sin(y_angle)  ; 0 1 0 ; -sin(y_angle) 0 cos(y_angle)];
    r_x = [1 0 0; 0 cos(x_angle) -sin(x_angle);0 sin(x_angle) cos(x_angle) ];
    r = r_x*r_y*r_z;
    for di = 1:Nx
        for dj = 1:Nx
           C2 = [ Range(di); Range(dj) ; 40e-3];
           XC2_m = zeros(size(XC2));
           for i_m=1:length(XC2)
            XC2_m(i_m,:) = ((r * XC2(i_m,:)') + C2)';
           end
           for ii=1:turn2
              B12 = zeros(size(P1));
              for j =1:N1
                 for i=1:N1
                     for jj=1:n_poly2
                        a = (j-1) + (i-1)*N1 + 1;
                        rage_Xc = (ii-1)*n_poly1+ 1 : (ii-1)*n_poly1 + n_poly1  + 1; 
                        if( is_on_poly(P1(:,a),XC1(rage_Xc,:),w) == 0  && is_inside_poly(P1(:,a),XC1,n_poly1) == 1 ) %not on the poly
                            B12(:,a) = B12(:,a) + integral(@(t)intgral_func2(t,XC2_m((ii-1) * n_poly2 +  jj ,:),XC2_m((ii-1) * n_poly2 +  jj + 1,:),P1(:,a),I),0,1,'ArrayValued',true);
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
                                M12(i_t,di,dj) = M12(i_t,di,dj) + B12(3,a);
                           end
                      end          
                  end
              end
            end
        end
    end
end
beep
M21 = M21 * dx1 * dx1 / I;
K = M21 / sqrt(L1*L2);
time = toc;