function [BCRB_CS,MSE_random]= Delay_estimation_random(radar,mu_alpha,sigma2_alpha,Tau,p_Tau,Psi,Phi,dim_MV,Nmc,SNR,snr_idx,N)

Nt = size(Psi,1);
tau_rand = zeros(1,Nmc);
tau = zeros(1,Nmc);
for mc_idx = 1:Nmc

    %---------------------------------------
    %     random time delay  uniform distribution
    tau(mc_idx) = radar.tau_min + (radar.tau_max-radar.tau_min) * rand(1);
    %random time delay
%     tau(mc_idx) = radar.tau_mean + radar.tau_sigma.*randn(1);
%     while tau(mc_idx)<radar.tau_min || tau(mc_idx)>radar.tau_max
%         tau(mc_idx) = radar.tau_mean + radar.tau_sigma.*randn(1);
%     end
    
%     %%%narrow the search
%     if SNR >= 25
%         narrow_center = floor(tau(mc_idx)/radar.dtau(snr_idx));
%         start = narrow_center - 10000;
%         stop = narrow_center + 10000;
%         if start < 0
%             start = 1;
%         end
%         if stop>length(Tau)
%             stop = length(Tau);
%         end
%         narrow_idx = start:stop;
%         Tau_temp = Tau;
%         p_Tau_temp = p_Tau;
%         Tau = Tau(narrow_idx);
%         p_Tau = p_Tau(narrow_idx);
%         Psi = transmit_waveform(radar,Nt,Tau);
%     end


    %standard Fisher information
    JF_CS(mc_idx) = JF_random(radar, mu_alpha, sigma2_alpha, Phi,tau(mc_idx),N);

%     [SNR,mc_idx]
%     %generate signal based on time delay
%     vchirp = exp(-1i*2*pi*(radar.B/2)*(radar.fasttime-tau(mc_idx))+1i*pi*radar.gamma*(radar.fasttime-tau(mc_idx)).^2).';
%     signal = (interp1(radar.pulse_time,radar.pulse_env,radar.fasttime-tau(mc_idx),'spline',0).').*vchirp;
%     %normalization
%     Etmp = radar.Ts*sum(abs(signal).^2);
%     signal = signal*sqrt(radar.Ep/Etmp)/sqrt(radar.Ep); 
%     %----------------------------------------------------
% 
%     %random noise
%     n = sqrt(radar.Pn) * sqrt(1/2) * (randn(Nt,1)+1i*randn(Nt,1));
%       %random amplitude
%     alpha = sqrt(sigma2_alpha) * sqrt(1/2)*(randn(1)+1i*randn(1)) + mu_alpha;
%     
%     %-------------------------------
%     %MMSE estimation, random mixing
%     %-------------------------------
% 
%     %generate a random sensing matrix
%     for i = 1:N
%         y_rand(:,i) = Phi(:,:,i)*(alpha*signal+n);
%     end
%     tau_rand(mc_idx) = MMSE_estimation_random(radar,alpha,mu_alpha,sigma2_alpha,Tau,p_Tau,Psi,Phi,y_rand,signal,dim_MV,N);
%  
%     if SNR >= 25
%         Tau = Tau_temp ;
%         p_Tau = p_Tau_temp;
%     end
end

JD_CS = mean(JF_CS);
JP = 0;
% JP = 1/(radar.tau_sigma)^2; %for Normal distribution
BCRB_CS = 1/(JD_CS+JP);

MSE_random = 0;
% BCRB_CS = 0;
% MSE_random = mean((tau_rand-tau).^2);

end

