function res = is_inside_poly(P,XC,n_poly)
res = 1;
for i=1:n_poly
   X1 = XC(i,:);
   X2 = XC(i+1,:); 
   if(X2(1) ~= X1(1) && X2(2) ~= X1(2))
       m = (X2(2) - X1(2)) / (X2(1) - X1(1));
       b = X1(2) - m * X1(1);
       xm = (X1(1) + X2(1)) / 2;
       ym = (X1(2) + X2(2)) / 2;
       if m >= 0 
           if xm < 0 
               if P(2)  >= ( m * P(1) + b )
                   res = 0;
                   break;
               end
           else
               if P(2)  <= ( m * P(1) + b )
                   res = 0;
                   break;
               end
           end
       else
           if xm < 0 
               if P(2)  <= ( m * P(1) + b )
                   res = 0;
                   break;
               end
           else
               if P(2)  >= ( m * P(1) + b )
                   res = 0;
                   break;
               end
           end
       end
   elseif X2(1) == X1(1)
       if X2(1) >= 0 && P(1) >= X2(1)
           res = 0;
           break;
       elseif X2(1) <= 0 && P(1) <= X2(1)
           res = 0;
           break;
       end
   elseif X2(2) == X1(2)
       if X2(2) >= 0 && P(2) >= X2(2)
           res = 0;
           break;
       elseif X2(2) <= 0 && P(2) <= X2(2)
           res = 0;
           break;
       end
   end
end
end

