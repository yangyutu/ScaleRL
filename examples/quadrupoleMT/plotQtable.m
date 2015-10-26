clear all
close all

x=linspace(0,1,20);
y=linspace(24000,16200,78);
[X Y] = meshgrid(x,y);

for i = 0 : 2
    filename = strcat('QTableFinal',num2str(i),'.dat');
    data{i+1} = load(filename);
    data{i+1} = data{i+1}';
    figure(i + 1)
    %imagesc(data{i+1});
    contourf(X,Y,data{i+1});
    colorbar;
    title('Q Value');
    xlabel('\psi_6')
    ylabel('Rg')
end

data{4} = load('actionMap.dat');
data{4}=data{4}';
data{4}(data{4} < 0) = NaN;
figure(5)
contourf(X,Y,data{4},'linestyle','none');
colorbar;