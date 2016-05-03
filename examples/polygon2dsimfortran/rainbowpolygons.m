clear all
close all
filestart=1;
fileend=100;


rgb=colormap('jet');
rgblen=size(rgb,1);
fileskip=2;
nfile=(fileend-filestart+1)/fileskip;

skip=2;
folderpath='I:\researchdata\polygon\assem\squ\3\';
%folderpath='N:\simulation data\quadrupole project\Bd dynamic\210p 0.7v\0.7v 210p starting from meshgrid stage 2\neworderwithlocalc6\op';
figure(1);
for i=filestart:fileskip:fileend
    filename=strcat(folderpath,'mc_order',num2str(i),'.txt');
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
    xlabel('psi');
    ylabel('Rg');
    hold on
    ylim([9.5 14])
    xlim([0 1])
    polish;
end
 tifname=strcat(folderpath,num2str(1),'.tif');

    print('-dtiff', tifname);
