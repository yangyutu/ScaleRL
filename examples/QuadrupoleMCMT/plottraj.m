clear all
close all
format long
folderpath='traj_cluster/';
lengthSet=[];
for i = 88
    data = load(strcat(folderpath,'controlop',num2str(i),'.dat'));
    figure(1)
    hold on
    plot(data(:,2),data(:,3));
    xlabel('\psi_6')
    ylabel('C_6')
    lengthSet = [lengthSet; data(end,1)];
    figure(11)
    plot(data(:,1),data(:,5),'linestyle','none','marker','o');
    xlabel('time, s')
    ylabel('action')    
    hold on
    figure(12)
    plot(data(:,1),data(:,2));
    xlabel('time, s')
    ylabel('\psi_6')    
    hold on
    
    
end

mean(lengthSet)

