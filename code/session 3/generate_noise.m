function [noise] = generate_noise(filename,signal,fs,SNR,h_noise,siglength)
[noise_raw{1}, noise_raw{2}] = audioread(filename);
noise_loudspeaker  = resample(noise_raw{1},fs,noise_raw{2});
noise_loudspeaker  = noise_loudspeaker(1:siglength*fs);
noise = fftfilt(h_noise,noise_loudspeaker);

noise_std = std(signal)/db2mag(SNR);
noise = noise_std./std(noise).*noise;
end
