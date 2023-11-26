clear
r = 1;
X = [3/2 , 0 , -3/2 , 0 , -3/2 , 0 , 3/2]*r;
Y = [-sqrt(3)/2 , -sqrt(3) , -sqrt(3)/2 , 0 , sqrt(3)/2 , sqrt(3) , sqrt(3)/2]*r;
Z = [ 3 , 0 , 0 , 0 , 0 , 0 , 0];
subplot(2,1,1)
scatter3(0,0,0)
hold on
h = scatterbar3(X,Y,Z,r/2);
h = squeeze(h);
axis equal

x_hex = [ r , r/2 , -r/2 , -r , -r/2 , r/2 , r];
y_hex = [ 0 , sqrt(3)/2 , sqrt(3)/2 , 0 , -sqrt(3)/2 , -sqrt(3)/2 , 0]*r;
for i =1:7
    plot(x_hex + X(i) , y_hex + Y(i) , 'black','linewidth',3)
end