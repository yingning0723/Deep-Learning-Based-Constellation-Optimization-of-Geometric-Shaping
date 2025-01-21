clear
figure
markersize =6 ;
linewidth = 0.9;
EbN0_db_qam = 5.5:0.2:7.7;
load("QAM64_without_gpu.mat")
ber_qam_withoutgpu = ber_qam;
load("ber_QAM64_gpu1.mat")
ber_qam_gpu1 = ber_qam;
load("ber_QAM64_gpu2.mat")
ber_qam_gpu2 = ber_qam2;
semilogy(EbN0_db_qam,ber_qam_withoutgpu,"Marker", 'o',"MarkerSize",markersize,"LineStyle","-","Color",[0.8 0.1 0.1],"LineWidth",linewidth);hold on
semilogy(EbN0_db_AE,ber_qam_gpu1,"Marker", 'p',"Markersize",markersize,"LineStyle","-","Color",[0.5 0.4 0.1],"LineWidth",linewidth);hold on
semilogy(EbN0_db_ATSC,ber_qam_gpu2,"Marker", 'square',"MarkerSize",markersize,"LineStyle","-","Color",[0.1 0.1 0.8],"LineWidth",linewidth);hold on
legend("BER-withoutGPU",'BER-gpuDecoder','BER-comm.gpu.LDPCDecoder',"Location","southwest")
xlabel("Eb/N0");ylabel("BER");
grid on
set(gca,'FontSize',14,'FontName','Times New Roman');
set(gca,'ytick',[1e-4,1e-3,1e-2,1e-1,1],'ygrid','on','gridlinestyle','-','Gridalpha',0.1)
ylim([1e-4,1])

hold off