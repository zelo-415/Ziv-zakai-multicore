clear;
clc;

CF = 10;%%压缩率
type = 0; %%Gaussian = 1; Uniform = 0
MMSE_switch = 0;
BCRB_switch = 0;
ZZB_switch = 1;
Phiopt_switch = 0;
%---------------------------------------
%%Fudemental Simulation Parameter
radar.c = 3e8;          %speed of light
radar.B = 5e6;          %Radar waveform bandwidth
radar.gamma = 6.25e11;  %Waveform chirp rate
%-----------------------------------------------------
%for compressive sensing
radar.B = radar.B*(CF/5);
radar.gamma = radar.gamma*(CF/5);
%-----------------------------------------------------
radar.waveform_type = 'chirp'; %radar waveform type
radar.Ts = 0.5/radar.B; %Underlying simulation sampling interval
radar.Fs = 1/radar.Ts; %Underlying simulation sample rate
radar.tau_min = 0; %minimum delay for interest

radar.tau_max = (40*(CF/5))/radar.B;

%Gaussian distribution of time delay
radar.tau_mean = mean([radar.tau_min,radar.tau_max]); %mean value of time delay
radar.tau_sigma = 5/(radar.B/(CF/5));
%---------------------------------------------------
radar.oversample = 4;
radar.waveform_filter_cut = 1/radar.oversample;
radar.waveform_fileter_order = 2*radar.oversample;
radar.waveform_filter = fir1(radar.waveform_fileter_order,2*radar.waveform_filter_cut);
%------------------------------------------------------
NTau_min = 100;
NTau_max = 1000;
% NTau_max = 10000;
% max_iter = 10000;
%------------------------------------------------------
%Number of Monte Carlo simulations
Nmc = 1000;


% %------------------------------------------------------

%SNR(dB)
% SNR = [-20:2:30,40,50];
% SNR = [40,50];
% SNR = [-20:2:10,24,30,40,50];
% SNR = 40;
SNR = 20;

%------------------------------------------------------
MSE_random = zeros(1,length(SNR));
for snr_idx = 1:length(SNR)
    
    SNR(snr_idx)
   
    %------------------------------------
    %Sensing matrix
    %-----------------------------------
    radar.SNR = 10.^(SNR(snr_idx)/10);
    
    mu_alpha = 0;
    sigma2_alpha = 1;
    
    radar.Ep = sigma2_alpha+mu_alpha^2;
    radar.N0 = 2*radar.Ep/radar.SNR;
    radar.Pn = (radar.N0/2)*radar.Fs;
    
    %Compute waveform
    radar = make_waveform(radar);
    
    Nt = ceil((radar.tau_max-radar.tau_min)/radar.Ts)+radar.Lp;
    radar.fasttime = radar.tau_min + (0:(Nt-1))*radar.Ts;
    
    if ZZB_switch == 1
        if CF~=10
            if SNR(snr_idx) < 20
                NTau_max = 200;
            elseif (SNR(snr_idx)>=20) && (SNR(snr_idx)<30)
                NTau_max = 2000;
            elseif (SNR(snr_idx)>=30) && (SNR(snr_idx)<40)
                NTau_max = 5000;
            elseif (SNR(snr_idx)>=40) && (SNR(snr_idx)<=50)
                NTau_max = 20000;
            else
                NTau_max = 100000;
            end
        elseif CF == 10
            if SNR(snr_idx) < 20
                NTau_max = 200;
            elseif (SNR(snr_idx)>=20) && (SNR(snr_idx)<25)
                NTau_max = 2000;
            elseif (SNR(snr_idx)>=25) && (SNR(snr_idx)<40)
                NTau_max = 10000;
            elseif (SNR(snr_idx)==40)
                NTau_max = 50000;
            elseif (SNR(snr_idx)==50)
                NTau_max = 200000;    
            end
        end
    end
    NTau_max = NTau_max/4;
    
    %here, Tau is for sensing matrix optimization
    Tau = linspace(radar.tau_min,radar.tau_max,NTau_max);
    %----------------------------------------------------------
    %generate transmit waveform for sensing matrix optimization
    Psi = transmit_waveform(radar,Nt,Tau);
    %dimension of measurement vector
    dim_MV = floor(Nt/(CF*radar.Fs/radar.B));
    %generate a random sensing matrix
%     Phi = random_sensing_matrix(Nt,dim_MV,radar.waveform_filter);
        K = 9;
%    for K = 1:10
%     for k = 1:K
%         Phi(:,:,k) = random_sensing_matrix(Nt,dim_MV,radar.waveform_filter);
%     end
%     save('Phi.mat','Phi')
    load('Phi.mat','Phi')
%     load('opt_Phi.mat','Phi');
%     load('opt_Phi_1.mat')
%     Phi = Phi_opt;
    %     %sensing matrix optimization
    if Phiopt_switch == 1
        if type == 0
            prior.pdf_type = 'uniform';
            NTau = length(Tau);
            p_Tau = 1/NTau * ones(1,NTau);
        elseif type == 1
            prior.pdf_type = 'Gaussian';
            p_Tau = normpdf(Tau,radar.tau_mean,radar.tau_sigma);
            p_Tau = p_Tau/sum1(p_Tau);
        end
        Phi_opt = opt_sensing_matrix_random(radar,p_Tau,Psi,Nt,dim_MV,CF,SNR(snr_idx),sigma2_alpha,mu_alpha);
    end
    %----------------------------------------------------
    
    if MMSE_switch == 1
        %---------------------------------------------
        %Cramer-Rao Bound
        %---------------------------------------------
        %classical CRB in frequency domain
        CRB(snr_idx) = CRB_delay_random(radar);
        
        radar.dtau(snr_idx) = sqrt(12*CRB(snr_idx)/10);
        %------------------------------------------------
        
        
        dtau = radar.dtau(snr_idx);
        Tau = radar.tau_min:dtau:radar.tau_max;
        NTau = length(Tau);
        if NTau <= NTau_min
            Tau = linspace(radar.tau_min,radar.tau_max,NTau_min);
            NTau = length(Tau);
        end
        
        if type == 0
            prior.pdf_type = 'uniform';
            NTau = length(Tau);
            p_Tau = 1/NTau * ones(1,NTau);
        elseif type == 1
            prior.pdf_type = 'Gaussian';
            p_Tau = normpdf(Tau,radar.tau_mean,radar.tau_sigma);
            p_Tau = p_Tau/sum1(p_Tau);
        end
        %         if SNR(snr_idx)<=25
        Psi = transmit_waveform(radar,Nt,Tau);
        %         end
    end
    
    
%     ------------------------------------------------
%     Phi = Phi_opt; 
%     -----------------------------------------------
  
    if BCRB_switch == 1
        [BCRB_CS(snr_idx),~]...
            =Delay_estimation_random(radar,mu_alpha,sigma2_alpha,Tau,0,Psi,Phi,dim_MV,Nmc,SNR(snr_idx),snr_idx,K);
    end
    
    if MMSE_switch == 1 
        [BCRB_CS(snr_idx),MSE_random(snr_idx)]...
            =Delay_estimation_random(radar,mu_alpha,sigma2_alpha,Tau,p_Tau,Psi,Phi,dim_MV,Nmc,SNR(snr_idx),snr_idx,K);
        save('MSE_zero_mean_exp_CF10_N8_uniform.mat','MSE_random','SNR')
    end
    
    
    if ZZB_switch == 1
            %---------------------------------------------
            %Ziv-Zakai Bound for Uniform Distribution
            %---------------------------------------------
            
            T = radar.tau_max;
            dT = radar.tau_max/NTau_max;
            Cnn = radar.Pn*eye(Nt);
            Cyy = radar.Pn*eye(dim_MV);
            iCyy = inv(Cyy);
            I = eye(dim_MV);
            Pe_total = 0;
            for tau_h = 127:NTau_max-1
                Pe_sum = 0;
                for tau_a = 1:NTau_max-tau_h
                    Psi_tau_a = Psi(:,tau_a)/sqrt(radar.Ep);
                    Psi_tau_ah = Psi(:,tau_a+tau_h)/sqrt(radar.Ep);
                    
                    C_psi_a = sigma2_alpha * Psi_tau_a*Psi_tau_a';
                    C_psi_ah = sigma2_alpha * Psi_tau_ah*Psi_tau_ah';
                    
                    for i = 1 : K
                        Cyy_tau_a(:,:,i) =  Phi(:,:,i) * C_psi_a * Phi(:,:,i)' + Cyy;
                        Cyy_tau_ah(:,:,i) =  Phi(:,:,i) * C_psi_ah * Phi(:,:,i)' + Cyy;
                        c(i) = log( abs( det( Cyy_tau_ah(:,:,i) ) / det( Cyy_tau_a(:,:,i) ) ) );
                    end
              
                    %
                    par = sum(c);
                   
                    
                    %%%第一类错误概率
                    
                    for i = 1 : K  
                        temp(:,:,i) = sqrtm( Cyy_tau_a(:,:,i) ) / Cyy_tau_ah(:,:,i) * sqrtm( Cyy_tau_a(:,:,i) ) - I;                      
                    end
                    temp0 = sum( temp , 3);
                  
                    if K > dim_MV/2
                        N = dim_MV/2;
                    else
                        N = K;
                    end
                    DD_a = sort(real(eigs(temp0,2*N)),'descend');
                    
                    eig_pos = sum(DD_a >= 0);
                    eig_neg = sum(DD_a <= 0);
%                
%                     u = 1./DD_a(1:N).';
%                     v = -1./DD_a(N+1:2*N).';
                    
                    u = 1./DD_a(1:eig_pos).';
                    v = -1./DD_a(eig_pos+1:2*N).';
                    
%                     c = ones(1,N);
%                     d = ones(1,N);
                    c = ones(1,eig_pos);
                    d = ones(1,eig_neg);
                    
%                     for j = 1:N                        
%                         c_j = u - u(j);
%                         d_j = v - v(j);                     
%                         c_j(c_j==0) = 1;
%                         d_j(d_j==0) = 1;                        
%                         c(j) = 1/prod(c_j);
%                         d(j) = 1/prod(d_j);             
%                     end
                    for j = 1:eig_pos
                       c_j = u - u(j);
                       c_j(c_j == 0) = 1;
                       c(j) = 1/prod(c_j);
                    end
                    for j = 1:eig_neg
                        d_j = v - v(j);
                        d_j(d_j == 0) = 1;
                        d(j) = 1/prod(d_j);
                    end
                    
%                     P1 = zeros(1,N);
%                     P2 = zeros(1,N);
                    P1 = zeros(1,eig_pos);
                    P2 = zeros(1,eig_neg);
                    
%                     for j = 1:N                   
%                         sum1 =  sum( d./ (u(j) + v) / u(j) );
%                         sum2 =  sum( c./ (u + v(j)) / v(j) );
%     
%                         P1(j) = c(j) * sum1;
%                         P2(j) = d(j) * sum2;
%                     end
                    for j = 1:eig_pos
                       sum1 = sum(d./ (u(j) + v) / u(j));
                       P1(j) = c(j) * sum1;
                    end
                    for j = 1:eig_neg
                       sum2 = sum(c./ (u + v(j)) / v(j));
                       P2(j) = d(j) * sum2;
                    end
                    K1 = prod(u);
                    K2 = prod(v);
                    
                    z = -par;
                    P1_sum = sum( P1 .* exp( -u * z ) );
                    P2_sum = sum( P2 .* exp( v * z ) );
                    
                    if z >= 0
                        Pe = 1 - K1 * K2 * P1_sum;
                    else
                        Pe = K1 * K2 * P2_sum;
                    end
 
 

%                     wa = real(DD_a/2);  
%                     wa = wa.'; 
%                     nc_a = zeros(1,2*N);
%                     dof = 2*ones(1,length(wa));
%                         
%                     ma = 0;
%                     
%                     s = 0;
%                     x = -par;
%                     
%                     Pe = gx2cdf(x,wa,dof,nc_a,ma,s,'AbsTol',0,'RelTol',1e-4);
                    %%%%第二类错误概率省略，因为两类错误对称
                    if Pe < 0 || Pe > 1 
                        continue;
                    end
                    Pe_sum = Pe_sum + Pe * dT;%%integration
%                     [SNR(snr_idx),tau_h,tau_a]
%                     tau_a
                end
                %%The second integration with dh
                    
                load laughter;
                sound(y,Fs);
                Pe_total = Pe_total + tau_h * dT * Pe_sum * dT;
                [SNR(snr_idx),tau_h]
                Pe_snr(snr_idx) = Pe_total/T
                if tau_h >= 2
                    if tau_h * dT * Pe_sum * dT / Pe_total < 1e-5
                        break;
                    end
                end
            end
            Pe_snr(snr_idx) = Pe_total/T;
%             save('ZZB_zero_mean_exp_CF10_N9_uniform.mat','Pe_snr','SNR')
    end

end
% sqrt(BCRB_CS)*1e6;
% temp = Pe_snr;
% SNR_temp = SNR
% save('BCRB_zero_mean_exp_CF10_N1_uniform.mat','BCRB_CS','SNR')