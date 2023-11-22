function Phi_opt = opt_sensing_matrix_random(radar,p_Tau,Psi,Nt,dim_MV,CF,SNR,sigma2_alpha,mu_alpha)
load('Phi.mat')
% Phi_init = random_sensing_matrix(Nt,dim_MV,radar.waveform_filter); %%%initial Phi is Gaussian random CS kernel
Phi_init = Phi(:,:,1);

Psi_bar = zeros(Nt,length(p_Tau));
Psi = Psi/sqrt(radar.Ep);
for k = 1:length(p_Tau)
    Psi_bar(:,k) = p_Tau(k).* Psi(:,k);
end
Psi_bar = sum(Psi_bar,2);

%------------------------------------------------
% gradient-based process
%-----------------------------------------------
% threshold = 1e-6;
Phi = Phi_init;
% gam = 0.2;
% if SNR <= 10
%     gam = 0.1;
% elseif SNR <=16
%     gam = 2e-3;
% elseif SNR <= 20
%     gam = 1e-4;
% if SNR <=30
% %     gam = 1e-6;%2e-5
% % gam = 2e-7;
%     gam = 0.01;
% elseif SNR<=50
%     gam = 1e-4;
% else
%     gam = 1e-5;
% end
% gam = -1e-6;

gam = 1e-4;
iter_num = 200;

% if mu_alpha ~= 0
%     for idx = 1:iter_num
%         [SNR,idx]
%         
%         expo = zeros(1,length(p_Tau));
%         q = zeros(1,length(p_Tau));
%         detCyy=zeros(1,length(p_Tau));
%         CyyNyq = zeros(dim_MV,Nt,length(p_Tau));
%         delta_Psi = zeros(length(Psi_bar),length(p_Tau));
%         Cyyk = zeros(dim_MV,dim_MV,length(p_Tau));
%        
%         I = eye(Nt);
%         Cnn = radar.Pn * Phi;
%         for k = 1:length(p_Tau)
%             CyyNyq(:,:,k) = (sigma2_alpha * Phi * Psi(:,k) * Psi(:,k)'+ Cnn);
%             Cyyk = reshape(CyyNyq(:,:,k),dim_MV,Nt) * Phi';
%             
%             delta_Psi(:,k) = Psi_bar-Psi(:,k);
%             q(k) = real(delta_Psi(:,k)'*Phi'/Cyyk*Phi*delta_Psi(:,k));
%             expo(k) = -mu_alpha^2*q(k);
%             detCyy(k) = det(Cyyk);
%         end
%         expo = expo - max(expo);
%         denominator = p_Tau.*exp(expo)./detCyy;
%         numerator = zeros(size(Phi,1),size(Phi,2));
%         term2 = zeros(size(Phi,1),size(Phi,2));
%         test = zeros(size(Phi,1),size(Phi,2));
%         for k = 1:length(p_Tau)
%             CyyNyq = Phi *(sigma2_alpha*Psi(:,k)*Psi(:,k)'+Cnn);
%             CyyNyq2 = reshape(CyyNyq(:,:,k),dim_MV,Nt);
%             Cyyk =  CyyNyq2 * Phi';
%            
%             fastter = Cyyk\CyyNyq2;
%             fastter2 =  Cyyk\Phi*delta_Psi(:,k)*delta_Psi(:,k)';
%             q_grad = 2*Cyyk\Phi*(Psi_bar-Psi(:,k))*(Psi_bar-Psi(:,k))'*(I-Phi'*fastter);
%             q_grad = fastter2 - fastter2 *Phi'*fastter;
%             Uk = fastter+mu_alpha^2*q_grad;
%             numerator = numerator + p_Tau(k)*exp(expo(k))/detCyy(k)*Uk;%(fastter+mu_alpha^2*q_grad);
%             term2 = term2 + p_Tau(k) *Uk;% (fastter+mu_alpha^2*q_grad);
%             test =  test + p_Tau(k)*q_grad;
%         end
%         grad = numerator/sum(denominator)-term2;
%         gam = 10^(floor(log10(abs(max(Phi(:)))/abs(max(grad(:))))));
%         
%         Phi = Phi + gam*grad;
%         
%         [~, ~, V] = svd(Phi);
%         Phi = V(1:dim_MV,:);
%         
%     end
    
% else
    for idx = 1:iter_num
        [SNR,idx]
        
        detCyy=zeros(1,length(p_Tau));
        tmp = zeros(dim_MV,Nt,length(p_Tau));
        tmp2 = zeros(dim_MV,Nt,length(p_Tau));
        Cnn = radar.Pn * eye(Nt);
        Cyy = radar.Pn * eye(dim_MV);
        for k = 1:length(p_Tau)
            CyyNyq = sigma2_alpha * Psi(:,k) * Psi(:,k)';
            CyyNyq =  CyyNyq + Cnn;
            Cyyk = sigma2_alpha * Phi * Psi(:,k) * Psi(:,k)' * Phi' + Cyy;
            detCyy(k) = det(Cyyk);
            tmp(:,:,k) = p_Tau(k) * Cyyk \ Phi * CyyNyq;
            tmp2(:,:,k) = tmp(:,:,k) / detCyy(k);
        end
        
        grad = sum(tmp2,3)/sum(p_Tau./detCyy) - sum(tmp,3);
        % gam = 10^(floor(log10(abs(max(Phi(:)))/abs(max(grad(:))))));
        
        Phi = Phi + gam*grad;
        
        [~, ~, V] = svd(Phi);
        Phi = V(1:dim_MV,:);
        
    end
% end

% [U, S, V] = svd(Phi);
% Phi = V(1:dim_MV,:);
Phi_opt = Phi;

save('opt_Phi_5.mat','Phi_opt')


