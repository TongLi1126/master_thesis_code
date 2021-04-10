function spl = psd2spl(signal,nfft)
%% calculate psd, spl(sound pressure of signal)
psdy = abs(signal).^2/(nfft*nfft); %%%%%%% got a problem here
%psdy = abs(signal).^2/(nfft*fs);
psdy(2:end-1) = 2*psdy(2:end-1); 
%plot(f,10*log10(psdy))
spl = 10*log10(psdy) +90;
% spl = 10*log10(signal) +104;
end

