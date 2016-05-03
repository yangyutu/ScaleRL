clear all
close all
filestart=1;
fileend=30;


rgb=colormap('jet');
rgblen=size(rgb,1);
fileskip=1;
nfile=(fileend-filestart+1)/fileskip;

skip=2;
figure(1);
filename='mc_order1.txt';
data=load(filename);

datalen=size(data,1);
psi=[];
Rg=[];
 for j=1:datalen
    if(mod(j,skip) ==1)
        psi=[psi;data(j,2)];
        Rg=[ Rg;data(j,3)];
    end
end
color=min(floor(rgblen/(fileend-filestart+1)*(i-filestart+1))+1,rgblen);
plot( psi,Rg,'color',rgb(color,:),'linewidth',3);
xlabel('S4');
ylabel('Rg');
hold on
ylim([7.2 9])
xlim([0 1])
 tifname=strcat(folderpath,num2str(1),'.tif');

    print('-dtiff', tifname);
