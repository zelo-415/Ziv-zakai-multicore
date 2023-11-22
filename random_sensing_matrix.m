%generate a random sensing matrix

function Phi = random_sensing_matrix(len_Psi,dim_MV,waveform_filter)

len_noise = dim_MV*len_Psi-length(waveform_filter)+1;

random_noise = randn(1,len_noise) + 1i*randn(1,len_noise);

%Filter the random noise
random_noise_cut = conv(random_noise, waveform_filter);

[U, S, V] = svd( reshape(random_noise_cut,len_Psi,length(random_noise_cut)/len_Psi).');

Phi = V(1:dim_MV,:); %Phi*Phi;=I




