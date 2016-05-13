clear
clc
close all
file = 'xyz_0.dat';
data = load(file);
pnum = max(data(:,1))+1;
data1 = data(1:2:end,:);
data2 = data(2:2:end,:);
dist = 1435*sqrt((data1(:,2)-data2(:,2)).^2+(data1(:,3)-data2(:,3)).^2)-2*1435;
meandist = mean(dist);
histhist = -hist(dist,100)/1000;
DMAX = max(dist);
DMIN = min(dist);
x = linspace(DMIN,DMAX,100);
plot(x,histhist)
figure(1)
% axis([0 5 -0.5 0])

% for i = 200:500
%     plot(data((i-1)*pnum+1:i*pnum,2),data((i-1)*pnum+1:i*pnum,3),'ko','markersize',6,'MarkerFaceColor','k')
%     axis([-30 30 -30 30])
%     title(strcat('t = ',num2str(i)))
%     if i == 1
%         xx = input('wait');
%     end
%     drawnow
%     mov(i) = getframe(gcf);
%     pause(0.1)
% end
% movie2avi(mov,'movie.avi','fps',5)
% range = data(601:900,:);
% Result = [];
% dist = ones(300,300);
% for i = 1:300
%     for j = 1:300
%         if (i~=j)
%             dist(i,j) = sqrt((data(i,2)-data(j,2))^2 + (data(i,3)-data(j,3))^2);
%              if dist(i,j) <= 2
%                 Result = [Result, [i;j]];
%              end
%         end
%     end
% end
% a = min(min(dist));