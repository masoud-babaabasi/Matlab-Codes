function [XC] = return_poly_spiral_corners(n_poly,turn,input_diameter,w,s)
R = input_diameter / 2;
XC = zeros( n_poly * turn + 1 , 3);
for i=1:turn
   for j=0:n_poly-1
       theta = 2 * pi / n_poly * j ;
       if j == 0
           XC( (i-1) * n_poly +  j + 1 , 1 ) = R * cos(theta) + (w+s)/sin(pi*(n_poly-2)/2/n_poly)/2;
           XC( (i-1) * n_poly +  j + 1 , 2 ) = R * sin(theta) - (w+s)/cos(pi*(n_poly-2)/2/n_poly)/2;
       else
           XC( (i-1) * n_poly +  j + 1 , 1 ) = R * cos(theta);
           XC( (i-1) * n_poly +  j + 1, 2 ) = R * sin(theta);
       end
   end
   R = R - ( w + s ) /cos(pi/n_poly);
end
theta = 0;
XC( end , 1 ) = R * cos(theta) + (w+s)/sin(pi*(n_poly-2)/2/n_poly)/2;
XC( end , 2 ) = R * sin(theta) - (w+s)/cos(pi*(n_poly-2)/2/n_poly)/2;
theta = 2*pi / n_poly / 2;
XCtx = XC(:,1) * cos(theta) - XC(:,2) * sin(theta);
XCty = XC(:,1) * sin(theta) + XC(:,2) * cos(theta);
XC(:,1) = XCtx;
XC(:,2) = XCty;
end