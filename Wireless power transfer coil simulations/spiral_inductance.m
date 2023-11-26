clear 
Rin = 7.5e-3 ; %Meter
w = 1e-3;
s = 0.35e-3;
N = 14;
I = 1;
C = [ 0 ; 0 ; 0 ];

R = Rin + N *(w + s);
% dx = 0.008 * R;% dS = 0.8%
dx = 1e-3;
Range = -(R-w/2):dx:(R-w/2);
Nr = length(Range);
P = [repmat(Range,1,Nr) ; kron(Range,ones(1,Nr)) ; zeros(1,Nr*Nr)];
B = zeros(size(P));
%%
for i=1:Nr
   for j =1:Nr
      a = (j-1) + (i-1)*Nr + 1;
      rp = sqrt((P(1,a)-C(1))^2 + (P(2,a)-C(2))^2);
          if( rp == 0 )
              tp = 0;
          else
            tp = atan( ( P(2,a) - C(2) ) / ( P(1,a) - C(1) ));
          end
          if( ( P(1,a) - C(1) ) < 0 )
              tp = tp + pi;
          end
          if( ( P(1,a) - C(1) ) > 0 && ( P(2,a) - C(2) ) < 0)
              tp = tp + 2*pi;
          end
%           for np = 0:N-1
%             tx = tp + 2*pi*np;
%             if(np == N-1)
%                 if rp < (Rin-w/2+(w+s)/2/pi*tx)
%                 con = 1;
%                 else 
%                     con = 0;
%                 end
%             else
%                 if (rp > (Rin+w/2+(w+s)/2/pi*tx) || rp < (Rin-w/2+(w+s)/2/pi*tx)) || P(3,a) ~= C(3)
%                     con = 1;
%                 else 
%                     con = 0;
%                     break;
%                 end
%             end
%           end

            if rp < (Rin-w/2+(w+s)/2/pi*tp)
                con = 1;
                else 
                    con = 0;
            end
          if con == 1
            B(:,a) = integral(@(theta)intgral_func3(theta,Rin,w,s,C,P(:,a),I),0,N*2*pi,'ArrayValued',true);
%              B(:,a) = [0,0,1];
          end
   end
end
L = sum(B(3,:)) * dx * dx / I;
%%
quiver3(P(1,:),P(2,:),P(3,:),B(1,:),B(2,:),B(3,:),'B')
hold on 
t = 0:1e-3:N*2*pi;
coil = [ (Rin + t * ( w + s ) / 2 /pi ) .*cos(t) ; (Rin + t * ( w + s ) / 2 /pi ) .*sin(t) ; zeros(1,length(t))];
plot3(coil(1,:),coil(2,:),coil(3,:),'r')