data = load('mc_xyz1.txt');
Pnum = max(data(:,1));
figure(1)
axis([-30 30 -30 30])

for i = 1:size(data,1)/Pnum
    plot(data((i-1)*Pnum+1:i*Pnum,2),data((i-1)*Pnum+1:i*Pnum,3),'ko', 'markersize',6);
    axis([-30 30 -30 30])
    drawnow
    title (strcat('t = ',num2str(i)))
    if i == 1
        x= input('wait');
    end
end

 
 