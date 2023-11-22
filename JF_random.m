


function JF = JF_random(radar, mu_alpha,sigma2_alpha,Phi,tau,N)

vchirp = exp(-1i*2*pi*(radar.B/2)*(radar.fasttime-tau)+1i*pi*radar.gamma*(radar.fasttime-tau).^2).';
signal = (interp1(radar.pulse_time,radar.pulse_env,radar.fasttime-tau,'spline',0).').*vchirp;

%normalization
Etmp = radar.Ts*sum(abs(signal).^2);
signal = signal * sqrt(radar.Ep/Etmp);

signal = signal/sqrt(radar.Ep);

Cnn = radar.Pn *eye(length(signal));
%CRB based on analytical derivation
% d_vchirp = (-1i*2*pi*(radar.B/2)-1i*2*pi*radar.gamma*(radar.fasttime-tau)).'.*vchirp;
d_vchirp = (1i*2*pi*(radar.B/2)-1i*2*pi*radar.gamma*(radar.fasttime-tau)).'.*vchirp;
d_signal = (interp1(radar.pulse_time,radar.pulse_env,radar.fasttime-tau,'spline',0).').*d_vchirp;

d_signal = d_signal/sqrt(radar.Ep);


for i = 1:N
    C_yy(:,:,i)  = Phi(:,:,i)  * (sigma2_alpha*signal*signal'+Cnn)*Phi(:,:,i)';
    C_yy_inv(:,:,i)  = inv(C_yy(:,:,i) );
    d_C_yy(:,:,i)  = Phi(:,:,i) * (d_signal * signal'+signal * d_signal')*Phi(:,:,i)';
    J(i) = real(signal'*Phi(:,:,i)'*C_yy_inv(:,:,i)*d_C_yy(:,:,i)*C_yy_inv(:,:,i)*Phi(:,:,i)*d_signal);
end

JF = 2 * sigma2_alpha * sum(J);





% function JF = JF_random(radar, mu_alpha,sigma2_alpha,Phi,tau,N)
% 
% vchirp = exp(-1i*2*pi*(radar.B/2)*(radar.fasttime-tau)+1i*pi*radar.gamma*(radar.fasttime-tau).^2).';
% signal = (interp1(radar.pulse_time,radar.pulse_env,radar.fasttime-tau,'spline',0).').*vchirp;
% 
% %normalization
% Etmp = radar.Ts*sum(abs(signal).^2);
% signal = signal * sqrt(radar.Ep/Etmp);
% 
% signal = signal/sqrt(radar.Ep);
% 
% Cnn = radar.Pn *eye(length(signal));
% C_yy = Phi * (sigma2_alpha*signal*signal'+Cnn)*Phi';
% C_yy_inv = inv(C_yy);
% 
% %CRB based on analytical derivation
% % d_vchirp = (-1i*2*pi*(radar.B/2)-1i*2*pi*radar.gamma*(radar.fasttime-tau)).'.*vchirp;
% d_vchirp = (1i*2*pi*(radar.B/2)-1i*2*pi*radar.gamma*(radar.fasttime-tau)).'.*vchirp;
% d_signal = (interp1(radar.pulse_time,radar.pulse_env,radar.fasttime-tau,'spline',0).').*d_vchirp;
% 
% d_signal = d_signal/sqrt(radar.Ep);
% 
% 
% % dd_vchirp = ((-1i*2*pi*(radar.B/2)-1i*2*pi*radar.gamma*(radar.fasttime-tau)).^2+1i*2*pi*radar.gamma).'.*vchirp;
% % dd_signal = (interp1(radar.pulse_time,radar.pulse_env,radar.fasttime-tau,'spline',0).').*dd_vchirp;
% % 
% d_C_yy = sigma2_alpha * Phi * (d_signal * signal'+signal*d_signal')*Phi';
% % dd_C_yy = sigma2_alpha*Phi*(dd_signal*signal'+2*d_signal*d_signal'+signal*dd_signal')*Phi';
% 
% JF = 2*sigma2_alpha*real(signal'*Phi'*C_yy_inv*d_C_yy*C_yy_inv*Phi*d_signal)...
%     + 2*mu_alpha^2*real(d_signal'*Phi'*C_yy_inv*Phi*d_signal);%-signal'*Phi'*C_yy_inv*dd_C_yy*C_yy_inv*Phi*signal);

