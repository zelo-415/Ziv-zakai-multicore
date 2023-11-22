clc;clear all;
load('ZZB_zero_mean_exp_CF10_N5_uniform.mat','Pe_snr','SNR');
load('BCRB_zero_mean_exp_CF10_N5_uniform.mat');
s = 1e6;
semilogy(SNR,sqrt(Pe_snr)*s,'linewidth',2);
grid on;hold on;
semilogy(SNR,sqrt(BCRB_CS)*s,'linewidth',2);
% load('BCRB_zero_mean_exp_CF10_N5_uniform2.mat');
% semilogy(SNR,sqrt(BCRB_CS)*s);
% load('opt_Phi_1.mat','Phi_opt')
% Phi(:,:,1) = Phi_opt;
% load('opt_Phi_2.mat','Phi_opt')
% Phi(:,:,2) = Phi_opt;
% load('opt_Phi_3.mat','Phi_opt')
% Phi(:,:,3) = Phi_opt;
% load('opt_Phi_4.mat','Phi_opt')
% Phi(:,:,4) = Phi_opt;
% load('opt_Phi_5.mat','Phi_opt')
% Phi(:,:,5) = Phi_opt;
% save('opt_Phi.mat','Phi');