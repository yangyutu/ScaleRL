data = load('xyz_0.dat');
Pnum = max(data(:,1))+1;
figure(1)
axis([-30 30 -30 30])
x= input('wait');
for i = 1:size(data,1)/Pnum
    plot(data((i-1)*Pnum+1:i*Pnum,2),data((i-1)*Pnum+1:i*Pnum,3),'ko', 'markersize',6);
    axis([-30 30 -30 30])
    drawnow
    title (strcat('t = ',num2str(i)))
    pause(0.5)
end

 
 