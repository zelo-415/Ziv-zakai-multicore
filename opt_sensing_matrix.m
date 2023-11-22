function Phi_opt = opt_sensing_matrix(radar,p_Tau,Psi,Nt,dim_MV,CF,SNR)
Phi_init = random_sensing_matrix(Nt,dim_MV,radar.waveform_filter); %%%initial Phi is Gaussian random CS kernel

alpha = 1;
Psi_bar = zeros(Nt,length(p_Tau));
for k = 1:length(p_Tau)
    Psi_bar(:,k) = p_Tau(k).*Psi(:,k);
end
Psi_bar = sum(Psi_bar,2);

%------------------------------------------------
% gradient-based process
%-----------------------------------------------
% threshold = 1e-6;
Phi = Phi_init;

gam = 1e-4;
iter_num = 500;
% gam = 1e-5;
% iter_num = 2000;
grad = zeros(size(Phi,1),size(Phi,2),2);

for idx = 1:iter_num    
    idx
    if idx == 1
        grad = zeros(size(Phi,1),size(Phi,2));
    end
    temp = grad;
    expo = zeros(1,length(p_Tau));
    for k = 1:length(p_Tau)
        expo(k) = -real(abs(alpha)^2/radar.Pn*(Psi_bar-Psi(:,k))'*Phi'*Phi*(Psi_bar-Psi(:,k)));
    end
    expo = expo - max(expo);
    denominator = p_Tau.*exp(expo);
    numerator = zeros(size(Phi,1),size(Phi,2));
    for k = 1:length(p_Tau)
        numerator = numerator + denominator(k).*Phi*(Psi_bar-Psi(:,k))*(Psi_bar-Psi(:,k))';
    end
    grad = abs(alpha)^2/radar.Pn*numerator/sum(denominator);
    
%     Phi =(1-gam)*Phi + gam*grad;
    Phi = Phi + gam*grad;
    [~, ~, V] = svd(Phi);
    Phi = V(1:dim_MV,:);
    delta(idx) = norm(temp,2)-norm(grad,2);
%     [norm(temp,2),norm(grad,2)]
%     delta(idx)
%     norm(grad)
%     if idx>2
% %         if delta(idx) > delta(idx-1) && delta(idx-1) >=0
% %             gam = gam/5;
% %         end
% %         if  delta(idx) < delta(idx-1) && delta(idx-1) <=0
% %             gam = gam*2;
% %         end
% %         if delta(idx) * delta(idx-1) < 0
% %             gam = gam/10;
% %         end
%         if delta(idx) <=threshold && delta(idx) >=0 && delta(idx)<delta(idx-1) 
%             break;
%         end
%     end     
end
% [U, S, V] = svd(Phi);
% Phi = V(1:dim_MV,:);
Phi_opt = Phi;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

% save(strcat('1022_phi_opt_uniform_CF',num2str(CF),'_SNR_',num2str(SNR)),'Phi_opt')


