function [SD,BSD,MBSD] = distor_cal(clean_stft,enhance_stft,nfft,TX_Mask)
% SD
x_STFT_norm = clean_stft./(max(max(max(clean_stft))));
S_mvdr_mwfL_stft_norm = enhance_stft./(max(max(max(enhance_stft))));
power=abs(x_STFT_norm(:,:,1))./abs(S_mvdr_mwfL_stft_norm(:,:,1));
power(find((power<=1e-6) |(power>=1e5))) = 1;
SD = mean(sqrt(mean((10*log10(power)).^2)));
% BSD and MBSD
LN_clean = psd2spl(x_STFT_norm(:,:,1),nfft);
LN_enhance = psd2spl(S_mvdr_mwfL_stft_norm(:,:,1) ,nfft);
L_diff = abs(LN_clean-LN_enhance); L_diff(find(L_diff>100)) = 0;
BSD = mean(sum((L_diff).^2))/mean(sum((LN_clean).^2));
M = (L_diff-TX_Mask); M(find(M>150)) = 0;
M(find(M>0)) = 1;M(find(M<=0)) = 0;
MBSD = mean(sum(M.*L_diff));
end