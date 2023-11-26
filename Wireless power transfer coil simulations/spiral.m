clear
fileName = 'spiral.bmp';
size = 1000;
N_turn = 10;
Dout = 1000;
Din = 300;
W = 10;

%**********************************
pic = ones(size,size);
alpha = ((Dout / 2) - ( Din/2 + W ) ) / (2 * pi * N_turn);
theta0 = (Din/2 + W/2) / alpha;

for theta = theta0 + (0:1e-5:2*pi*N_turn)
    for i=-W/2:W/2; 
        r = alpha * theta + i ;
        x = floor( r * cos(theta));
        y = floor( r * sin(theta));
        pic(x + size/2, y + size /2) = 0;
    end
end
imwrite(pic, fileName);