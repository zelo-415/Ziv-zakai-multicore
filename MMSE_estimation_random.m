function tau_MMSE = MMSE_estimation_random(radar,alpha,mu_alpha,sigma2_alpha,Tau,p_Tau,Psi,Phi,y,signal,dim_MV,N)

[Nt, K] = size(Psi);
C_nn = radar.Pn*eye(Nt);
f_y_k = zeros(N,K);
f_Y_k = zeros(1,K);
%---------------------------------
%propbability of measurement given time delay
for k = 1:K
   for i = 1:N  
        C_yy_Phi_k(:,:,i) = Phi(:,:,i) * (sigma2_alpha * Psi(:,k) * Psi(:,k)'+ C_nn) * Phi(:,:,i)';
        f_y_k(i,k) = real(log(1/det(C_yy_Phi_k(:,:,i)))-(y(:,i)-mu_alpha*Phi(:,:,i)*Psi(:,k))'/(C_yy_Phi_k(:,:,i))*(y(:,i)-mu_alpha*Phi(:,:,i)*Psi(:,k)));
   end
end
f_Y_k = sum(f_y_k,1);
scale = max(f_Y_k);
f_Y_k = exp(f_Y_k-scale);

%--------------------------------------

%MMSE estimate

tau_MMSE = sum(p_Tau.*Tau.*f_Y_k)/sum(p_Tau.*f_Y_k);

end
