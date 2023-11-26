function res = is_on_poly(P,XC,w)
angle = atan( P(2) / P(1) );
if P(2) > 0 && P(1) < 0
    angle = angle + pi;
elseif P(2) < 0 && P(1) < 0 
    angle = angle + pi;
elseif P(2) < 0 && P(1) > 0 
    angle = angle + 2*pi;
end
res = 0;
for i=1:length(XC)-1
    X1 = XC(i,:);
    X2 = XC(i+1,:);
    a = X2(1) - X1(1);
    b = X2(2) - X1(2);
    angle1 = atan( X1(2) / X1(1) );
    if X1(2) > 0 && X1(1) < 0
        angle1 = angle1 + pi;
    elseif X1(2) < 0 && X1(1) < 0 
        angle1 = angle1 + pi;
    elseif X1(2) < 0 && X1(1) > 0 
        angle1 = angle1 + 2*pi;
    end
    
    angle2 = atan( X2(2) / X2(1) );
    if X2(2) > 0 && X2(1) < 0
        angle2 = angle2 + pi;
    elseif X2(2) < 0 && X2(1) < 0 
        angle2 = angle2 + pi;
    elseif X2(2) < 0 && X2(1) > 0 
        angle2 = angle2 + 2*pi;
    end
    distance = abs( a * (X1(2) - P(2)) - (X1(1) - P(1)) * b )/sqrt(a^2 + b^2);
    if angle1 < angle2 
        if((distance < (w / 2)) && (angle >= angle1 && angle <= angle2) )
            res = 1;
            break;
        end
    else
        if((distance < (w / 2)) && (angle >= angle1 || angle <= angle2) )
            res = 1;
            break;
        end
    end
end
end

