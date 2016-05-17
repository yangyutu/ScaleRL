clear
clc
close all
file = 'xyz_0.dat';
data = load(file);
pnum = max(data(:,1))+1;
dist = [];
data1 = data(1:2:end,:);
data2 = data(2:2:end,:);
ddist = 1435*sqrt((data1(:,2)-data2(:,2)).^2+(data1(:,3)-data2(:,3)).^2)-2*1435;
for i = 1:10000
    if ddist(i) < 500
        dist = [dist ddist(i)];
    end
end
% meandist = mean(dist);
% y = histogram(dist,200);

histhist = hist(dist,50)/9634;
U = -log(histhist);
for i = 1:size(U,2)
    if U(i) >10
        U(i) = 0;
    end
end
Umax = max(U);
for i = 1:size(U,2)
    if U(i) == 0
        U(i) = Umax;
    end
end
U = U - max(U);
DMAX = max(dist);
DMIN = min(dist);
x = linspace(DMIN,DMAX,50);
plot(x,U)
% figure(1)
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