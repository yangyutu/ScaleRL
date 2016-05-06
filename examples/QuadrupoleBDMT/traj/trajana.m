clear
clc
close all
file = 'xyz_0.dat';
data = load(file);
pnum = max(data(:,1))+1;
for i = 1:size(data,1)/pnum
    plot(data((i-1)*pnum+1:i*pnum,2),data((i-1)*pnum+1:i*pnum,3),'ko','markersize',6)
    axis([-30 30 -30 30])
    title(strcat('t = ',num2str(i)))
%     if i == 1
%         xx = input('wait');
%     end
%     drawnow
    pause(0.5)
end
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