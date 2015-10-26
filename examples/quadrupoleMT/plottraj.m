clear all
close all
format long
folderpath='traj/';
lengthSet=[];
for i = 0 : 1: 100
    data = load(strcat(folderpath,num2str(i),'.dat'));
    figure(1)
    hold on
    plot(data(:,2),data(:,3));
    lengthSet = [lengthSet; length(data(:,2))];
end

lengthSet2=[];
for i = 0 : 99
    data = load(strcat(folderpath,'quench',num2str(i),'.dat'));
    figure(2)
    hold on
    plot(data(:,2),data(:,3));
        lengthSet2 = [lengthSet2; length(data(:,2))];
end
mean(lengthSet)
max(lengthSet)
mean(lengthSet2)