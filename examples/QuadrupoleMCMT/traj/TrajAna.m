data = load('xyz_0.dat');
Pnum = max(data(:,1))+1;
figure(1)
axis([-30 30 -30 30])
dist = 5*ones(100,100);
for i = 1:size(data,1)/Pnum
    plot(data((i-1)*Pnum+1:i*Pnum,2),data((i-1)*Pnum+1:i*Pnum,3),'ko', 'markersize',6);
    for ii = 1:100
        for jj = 1:100
            if ii ~= jj
                dist(ii,jj) = sqrt((data((i-1)*Pnum+ii,2)-data((i-1)*Pnum+jj,2))^2+(data((i-1)*Pnum+ii,3)-data((i-1)*Pnum+jj,3))^2);
            end
        end
    end
    k(i) = min(min(dist));
    axis([-30 30 -30 30])
    drawnow
    title (strcat('t = ',num2str(i)))
    if i == 1
        x= input('wait');
    end
    pause(0.5)
end


 
 