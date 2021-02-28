%%   Implementation of the MWF
num_mics = 2;
Rnn_pre = cell(N_freqs,1);  Rnn_pre(:) = {1e-6*ones(num_mics,num_mics)};      % Noise Only (NO) corr. matrix. Initialize to small random values
Ryy_pre = cell(N_freqs,1);  Ryy_pre(:) = {1e-6*ones(num_mics,num_mics)};
Rnn = cell(N_freqs,1);  Rnn(:) = {1e-6*ones(num_mics,num_mics)};      % Noise Only (NO) corr. matrix. Initialize to small random values
Ryy = cell(N_freqs,1);  Ryy(:) = {1e-6*ones(num_mics,num_mics)};      % Speech + Noise (SPN) corr. matrix. Initialize to small random values
lambda = 0.995;                                                       % Forgetting factors for correlation matrices - can change
SPP_thr = 0.95;                                                       % Threshold for SPP - can change

% For MWF filter for left ear
S_mvdr_mwfL_stft = zeros(N_freqs,N_frames,num_mics);
X_mvdr_mwfL_stft = zeros(N_freqs,N_frames,num_mics);
N_mvdr_mwfL_stft = zeros(N_freqs,N_frames,num_mics);
W_mvdr_mwfL = (1/num_mics)*ones(num_mics,N_freqs,N_frames);
I = eye(num_mics);
% STFT Processing
% Looping through each time frame and each frequency bin
tic
for l=1:N_frames % Time index
    
    for k = 1:N_freqs % Freq index
        
        
        % Create a vector of mic. signals
        Y_kl = squeeze(y_STFT(k,l,1:num_mics));  % M x 1 noisy mic sigs for this freq and time frame
        X_kl = squeeze(x_STFT(k,l,1:num_mics));
        N_kl = squeeze(n_STFT(k,l,1:num_mics));
        
        % Calculate alpha and m
        alpha = 0;
        m = 0.5;
        
        % ## Update the correlation matrices using the forgetting factor.
        % Threshold the SPP in order to distinguish between periods of speech and non-speech activity
%         if(SPP(k,l)>SPP_thr)
%             Ryy_pre{k} = lambda^2*Ryy_pre{k}+(1-lambda^2).*Y_kl*Y_kl';
%             % speech + noise
%         else
%             % noise only
%             Rnn_pre{k} = lambda^2*Rnn_pre{k}+(1-lambda^2).*Y_kl*Y_kl'; % N*N' is not available in real situation 
%         end
% %         Ryy{k} = (1-alpha)*Ryy_pre{k}+alpha*I;
% %         Rnn{k} = (1-alpha)*Rnn_pre{k}+alpha*(1-m)*I;
%           Ryy{k} = Ryy_pre{k};
%           Rnn{k} = Rnn_pre{k};

        if(SPP(k,l)>SPP_thr)
            Ryy{k} = lambda^2*Ryy{k}+(1-lambda^2).*Y_kl*Y_kl';
            % speech + noise
        else
            % noise only
            Rnn{k} = lambda^2*Rnn{k}+(1-lambda^2).*Y_kl*Y_kl'; % N*N' is not available in real situation 
        end
       % mu1: fix constant; mu2 spp model; mu3 psycho1;mu4 psycho2
       mu = mu_calc(SPP(k,l),TX_Mask(k,l),NMR(k,l),1);
       muxx(k,l) =  mu;
             
        % Computing the MWF filter using a generalised eigenvalue
        % decomposition of the correlation matrices.
      [V,d] = eig(Ryy{k},Rnn{k},'vector');
      [d,ind] = sort(d,'descend');
      V = V(:,ind);
      Qh = pinv(V);  
     sig_yy = V'*Ryy{k} *V;
sig_yy = [sig_yy(1,1) 0;0 sig_yy(2,2)];
sig_nn = V'*Rnn{k}*V;
sig_nn = [sig_nn(1,1) 0;0 sig_nn(2,2)];
Q= Qh';
sig_ss = (sig_yy(1,1) - sig_nn(1,1)) * Q(:,1)*Qh(1,:); 
W_update1 = pinv(Ryy{k})*sig_ss *[1;0];
weight_model = 3;power_model = 2;
W_update2 = weight_cal(Ryy{k}, Rnn{k},power_model,weight_model,mu);    
%       W = cal_w(Ryy{k},Rnn{k});
            % See DASP-CH4-37/40
            % d(1) corresponds to (d_noise/d_s&n)
            W_update = (1-1/d(1))*Qh(1,1).*V(:,1);% Final expression for filter.
%           W_update = (1-mu/(d(1)+mu-1))*Qh(1,1).*V(:,1);



            % After experiments, we found that Filter parameter is not likely
            % larger than 2 if it converges.
            if(max(abs(W_update)) < 2.2)  % This will skip some initial frames  
                W_mvdr_mwfL(:,k) =W_update1;  
                weight1(:,k,l) = W_update1;
                weight(:,k,l) = W_update;
            end
%            W_mvdr_mwfL1(:,k) =W_update1;
            % Filtering the noisy speech, the speech-only, and the noise-only.
        for m = 1:2
             S_mvdr_mwfL_stft(k,l,m) = W_mvdr_mwfL(:,k)'* Y_kl(1:num_mics);
%              X_mvdr_mwfL_stft(k,l,m) = W_mvdr_mwfL(:,k)'* X_kl(1:num_mics);
             N_mvdr_mwfL_stft(k,l,m) = W_mvdr_mwfL(:,k)'* N_kl(1:num_mics);
        end
        
        
    end % end freqs
end % end time frames
toc
% Observe processed STFTst

figure; subplot(2,1,1);
imagesc(time,f/1000,mag2db(abs(y_STFT(:,:,1))), [clow, chigh]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('microphne signal, 1st mic');
S_L_enhanced = S_mvdr_mwfL_stft(:,:,1);
subplot(2,1,2); imagesc(time,f/1000,mag2db(abs(S_mvdr_mwfL_stft(:,:,1))), [clow, chigh]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'),title('Enhanced Signal L - MWF');
% figure;
% imagesc(1:N_frames, f,muxx); colorbar; axis xy; set(gcf,'color','w');
% set(gca,'Fontsize',14), xlabel('Time Frames'), ylabel('Frequency (Hz)'), title('Masking threshold for ref mic');

% Apply the synthesis stage of the WOLA framework to obtain the time domain equivalents:

s_mwfL = WOLA_synthesis(S_mvdr_mwfL_stft,window,nfft,noverlap);% To complete (time-domain version of S_mvdr_mwfL_stft)
x_mwfL = WOLA_synthesis(X_mvdr_mwfL_stft,window,nfft,noverlap);% To complete (time-domain version of X_mvdr_mwfL_stft)
n_mwfL = WOLA_synthesis(N_mvdr_mwfL_stft,window,nfft,noverlap);% To complete (time-domain version of N_mvdr_mwfL_stft)


% PLOT SIGNALS
figure;
subplot(2,1,1);
plot(noisy_sig); title('signal before enhancement');
subplot(2,1,2);
plot(real(s_mwfL));title('signal after enhancement');
%  soundsc(s_mwfL,fs);
% soundsc(noisy_sig,fs_RIR);
% LISTEN TO SIGNALS!
%soundsc(s_mwfL,fs);
% soundsc(y_TD,fs);

% EVALUATION

%  SNR_in_L =snr(speech_non_rev(:,1),noise(:,1)) % Compute input SNR
% SNR_out_L =snr(s_mwfL(:,1),n_mwfL(:,1)) % Compute output SNR
% delta_SNR_L = SNR_out_L-SNR_in_L
%  SNR_in_R =snr(speech_non_rev(:,2),noise(:,2)) % Compute input SNR
% SNR_out_R =snr(s_mwfL(:,2),n_mwfL(:,2))
% delta_SNR_R = SNR_out_R-SNR_in_R
SNR_in_L =snr(speech(:,1),noise(:,1)); % Compute input SNR
SNR_out_L =snr(s_mwfL(:,1),n_mwfL(:,1)) ;% Compute output SNR
delta_SNR_L = SNR_out_L-SNR_in_L
% SNR_in_R =snr(speech_rev(:,2),noise_rev(:,2)) % Compute input SNR
% SNR_out_R =snr(s_mwfL(:,2),n_mwfL(:,2))
% delta_SNR_R = SNR_out_R-SNR_in_R
% spectral distance
power_L=abs(y_STFT(:,:,1))./abs(S_mvdr_mwfL_stft(:,:,1));
power_L = power_L(:,3:N_frames);
SD_L = sum(sqrt(sum((10*log10(power_L)).^2,1)/N_freqs),2)/N_frames