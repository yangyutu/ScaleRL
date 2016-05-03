% data = load('mc_xyz1.txt');
Pnum = 300;
for i = 1:size(data,1)/Pnum
    plot(data((i-1)*Pnum+1:i*Pnum,2),data((i-1)*Pnum+1:i*Pnum,3),'ko');
    axis([-30 30 -30 30])
    title (strcat('t = ',num2str(i)))
end