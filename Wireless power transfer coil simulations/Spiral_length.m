function [L] = Spiral_length(dout,w,s,n,num)
L = zeros(1,num);
for j =1:num
    m =  4 * n(j) ;

    for i=1:m
       if(i <= 3) 
          L(1,j) = L(1,j) + dout;
       else
           if(mod(i,2) == 0 )
               L(1,j) = L(1,j) + dout - (i-2)/2*(w + s);
           else L(1,j) = L(1,j) + dout - (i-3)/2*(w + s);
           end
       end
    end
end
end

