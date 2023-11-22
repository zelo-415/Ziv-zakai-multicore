load('BCRB_zero_mean_exp_CF10_N1_uniform.mat')
semilogy(SNR,sqrt(BCRB_CS)*1e6,'k--','Linewidth',1.5);
hold on; grid on;
load('ZZB_zero_mean_exp_CF10_N1_uniform.mat')
semilogy(SNR,sqrt(Pe_snr)*1e6,'k','Linewidth',1.5);

load('BCRB_zero_mean_exp_CF10_N2_uniform.mat')
semilogy(SNR,sqrt(BCRB_CS)*1e6,'g--','Linewidth',1.5);

load('ZZB_zero_mean_exp_CF10_N2_uniform.mat')
semilogy(SNR,sqrt(Pe_snr)*1e6,'g','Linewidth',1.5);

load('BCRB_zero_mean_exp_CF10_N4_uniform.mat')
semilogy(SNR,sqrt(BCRB_CS)*1e6,'b--','Linewidth',1.5);

load('ZZB_zero_mean_exp_CF10_N4_uniform.mat')
semilogy(SNR,sqrt(Pe_snr)*1e6,'b','Linewidth',1.5);

load('BCRB_zero_mean_exp_CF10_N8_uniform.mat')
semilogy(SNR,sqrt(BCRB_CS)*1e6,'r--','Linewidth',1.5);

load('ZZB_zero_mean_exp_CF10_N8_uniform.mat')
semilogy(SNR,sqrt(Pe_snr)*1e6,'r','Linewidth',1.5);

set(gca,'FontSize',14,'Fontname', 'Times New Roman');
xlabel('SNR (dB)','FontSize',14,'FontName','Times New Roman');
ylabel('RMSE (\mus)','FontSize',14,'FontName','Times New Roman');

legend('BCRB(K=1)','ZZB(K=1)','BCRB(K=2)','ZZB(K=2)','BCRB(K=4)','ZZB(K=4)','BCRB(K=8)','ZZB(K=8)','FontSize',14,'FontName','Times New Roman');


%%%%SNR = 50 dB, varying K
% load('BCRB_zero_mean_CF10_SNR50K1to8_uniform.mat')
% semilogy(sqrt(BCRB_CS)*1e6,'Linewidth',1.5);
% hold on; grid on;
% load('ZZB_zero_mean_CF10_SNR50K1to8_uniform.mat')
% semilogy(sqrt(Pe_snr)*1e6,'Linewidth',1.5);
% set(gca,'FontSize',14,'Fontname', 'Times New Roman');
% xlabel('K','FontSize',14,'FontName','Times New Roman');
% ylabel('RMSE (\mus)','FontSize',14,'FontName','Times New Roman');
% legend('BCRB','ZZB','FontSize',14,'FontName','Times New Roman');
