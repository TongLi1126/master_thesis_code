function [] = illustrate (RT_set,SNR_set,data,title)
figure;
for line = 1:length(SNR_set)+1
    if line ==1       
    plot(RT_set,data(:,line),'-s','MarkerSize',10, 'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
    else
   plot(RT_set,data(:,line),'-s');     
    end
    hold on;
end
legend('SNR =-5, TX','SNR =-5,MWF','SNR =0,MWF','SNR =5,MWF','SNR =10,MWF','SNR =15,MWF','SNR =20,MWF');
grid on;ylabel(title);xlabel('Reverberation time T_{60}')
% saveas(gcf,['/users/students/r0771119/Documents/master_thesis/master_thesis_code/result',title,'.png']);
end