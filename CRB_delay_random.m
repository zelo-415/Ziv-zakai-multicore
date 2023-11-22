
%CRB for time delay
function CRB = CRB_delay_random(radar)

vchirp = exp(-1i*2*pi*(radar.B/2)*(radar.fasttime)+1i*pi*radar.gamma*(radar.fasttime).^2).';
xt_tmp = radar.pulse_env.*vchirp(1:radar.Lp);
Xf = radar.Ts*abs(fftshift(fft(xt_tmp)));
Xf2 = Xf.^2;
df = 1/(radar.Lp*radar.Ts);
freq = -0.5*radar.Fs +(0:(radar.Lp-1))*df;
f_mean = df*sum(Xf2.*freq');
radar.B_rms = df*sum(Xf2.*(freq'-f_mean).^2)*(2*pi)^2;

%The CRB of signal delay
CRB = 1/(radar.SNR*radar.B_rms);

