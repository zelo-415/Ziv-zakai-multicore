
%for fixed amplitude
function Psi = transmit_waveform(radar, Nt, Tau)

NTau = length(Tau);
Psi = zeros(Nt,NTau);
%Cnn = zeros(Nt,Nt,NTau);

for idx = 1:NTau
    
    vchirp = exp(-1i*2*pi*(radar.B/2)*(radar.fasttime-Tau(idx))+ 1i*pi*radar.gamma*(radar.fasttime-Tau(idx)).^2).';
    Psi(:,idx) = (interp1(radar.pulse_time,radar.pulse_env,radar.fasttime-Tau(idx),'spline',0).').*vchirp;
    %normalization
    Etmp = radar.Ts*sum(abs(Psi(:,idx)).^2);
    Psi(:,idx) = Psi(:,idx)*sqrt(radar.Ep/Etmp);

    %Cnn(:,:,idx)=radar.Pn*eye(Nt);
end

