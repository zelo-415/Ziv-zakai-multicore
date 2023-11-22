function radar = make_waveform(radar)

radar.Tp = radar.B/radar.gamma; %radar waveform pulsewidth
radar.BT = radar.Tp*radar.B;%radar waveform time-bandwidth product

%Make the finely sampled radar pulse envelope with smooth transitions
L1 = ceil(radar.Tp/radar.Ts);
radar.pulse_env = conv(hanning(fix(0.1*L1)),ones(L1,1));
radar.Lp = length(radar.pulse_env);
M1 = 0.5*radar.Lp;
pulse_time = (-M1:M1)'*radar.Ts + 0.5 *radar.Tp;
radar.pulse_time = pulse_time(1:length(radar.pulse_env));

%Normalize the energy
Etmp = radar.Ts*sum(abs(radar.pulse_env).^2);
radar.pulse_env = radar.pulse_env*sqrt(radar.Ep/Etmp);

tmp_idx = find( (abs(radar.pulse_env)-0.5*max(abs(radar.pulse_env))) >0 );
radar.Lp_eff = fix(max(tmp_idx)-min(tmp_idx));
radar.waveform = radar.pulse_env.*exp(-1i*pi*radar.B*radar.pulse_time+1i*pi*radar.gamma*radar.pulse_time.^2).';

end

