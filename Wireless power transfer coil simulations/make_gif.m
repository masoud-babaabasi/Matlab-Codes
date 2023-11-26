clear 
load 'mutual3.mat' K
clear P1 P2
R = 25e-3 ; %Meter
w = 1e-3;
s = 0.35e-3;
I = 1;
% dx = 10e-3;
dx = R/3;
err = 10e-3;
dx = dx / R;
err = err / R;
w = w / R;
s = s / R;
R = 1;
D = -3:0.2:3;
Di = -3:0.05:3;

K0 = squeeze(K(4,D == 0,:))*100;
K0i = interp1(D,K0,Di);
fig = figure;
gif_delay = 0.01;

for idx = 1:length(K0i)
    subplot(2,1,1)
    marker_size = 10;
    plot(R,0,'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',marker_size)
    hold on
    plot(-R,0,'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',marker_size)

    xlim([-3 3])
    ylim([-1 3])
    xlabel('Y')
    ylabel('Z')
    
    plot(Di(idx),1.1,'o','MarkerFaceColor','blue','MarkerEdgeColor','blue','MarkerSize',marker_size)
    plot(Di(idx),2.5,'o','MarkerFaceColor','blue','MarkerEdgeColor','blue','MarkerSize',marker_size)
    hold off
    
    subplot(2,1,2)
    plot(Di(1:idx),K0i(1:idx),'LineWidth',3)
    xlim([D(1) D(end)])
    ylim([min(K0) max(K0)])
    hold on
    yline(0);
    set(gca,'FontSize',16,'FontWeight','Bold')
    hold off
    drawnow
    frame = getframe(fig);
    im{idx} = frame2im(frame);
end
filename = 'K_angle_0.gif'; % Specify the output file name
for idx = 1:length(K0i)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',1,'DelayTime',gif_delay);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',gif_delay);
    end
end