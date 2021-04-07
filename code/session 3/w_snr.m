function w_snr =  w_snr(x_stft,n_stft,bark)
    N_freqs = length(bark);
    I = [0;0.0103;0.0261;0.0419;0.0577;0.0577;0.0577;0.0577;0.0577;0.0577;
         0.0577;0.0577;0.0577;0.0577;0.0577;0.0577;0.0577;0.0577;0.0460;
         0.0343;0.0226;0.0110];
% I = ones(22,1)./22;
    N_frams = size(x_stft,2);
    i = 1; criti_band = [];
    sig_pow = zeros(21,N_frams);noise_pow = zeros(21,N_frams);snr =zeros(21,N_frams);
    for k = 2:N_freqs

        if(bark(k)>i)
            snr(i,:) = 10*log( sig_pow(i,:)./ noise_pow(i,:));
            i=i+1; 
            sig_pow(i,:) = 0; 
            noise_pow(i,:) = 0; 
            snr(i,:) = 0;
        end
         sig_pow(i,:) = sig_pow(i,:)+abs(x_stft(k,:)).^2; 
         noise_pow(i,:) = noise_pow(i,:)+abs(n_stft(k,:).^2);          
         criti_band = [criti_band;k i];    
    end
    w_snr = sum(I.*snr);
    w_snr(find(w_snr <=-100)) = [];
    w_snr = w_snr(~isnan(w_snr));
    w_snr=max(-20,w_snr); w_snr=min(20,w_snr);
    w_snr = mean(w_snr);

end